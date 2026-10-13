#!/usr/bin/env python3
"""Versioned SBDB POI selection/enrichment and Horizons SPK extraction.
Standard library only. Requests are sequential, cached, and resumable.
"""
import argparse
import base64
import datetime as dt
import hashlib
import json
import math
import sys
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request

SBDB_QUERY = 'https://ssd-api.jpl.nasa.gov/sbdb_query.api'
SBDB = 'https://ssd-api.jpl.nasa.gov/sbdb.api'
HORIZONS = 'https://ssd.jpl.nasa.gov/api/horizons.api'
FIELDS = 'spkid,pdes,full_name,name,prefix,class,diameter,diameter_sigma,albedo,H,GM,density,rot_per,pole,spec_B,spec_T,epoch,e,a,q,i,om,w,ma,tp,condition_code,orbit_id'
COMET_FIELDS = FIELDS + ',M1,K1,M2,K2,A1,A2,A3'
DEFAULT_MANIFEST = Path(__file__).with_name('selection.json')


def timestamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


class Client:
    def __init__(self, cache, offline=False, refresh=False, delay=0.2):
        self.cache, self.offline, self.refresh, self.delay = Path(cache), offline, refresh, delay
        self.last_request = 0.
        self.sources = {}

    def get(self, endpoint, parameters):
        url = endpoint + '?' + urllib.parse.urlencode(sorted(parameters.items()))
        key = hashlib.sha256(url.encode()).hexdigest()
        path = self.cache / (key + '.json')
        if path.exists() and not self.refresh:
            envelope = json.loads(path.read_text())
            if envelope.get('url') != url or envelope.get('schema') != 1:
                raise ValueError('Invalid cache envelope: ' + str(path))
        else:
            if self.offline:
                raise ValueError('Offline cache miss: ' + url)
            for attempt in range(4):
                time.sleep(max(0, self.delay - (time.monotonic() - self.last_request)))
                self.last_request = time.monotonic()
                request = urllib.request.Request(url, headers={'User-Agent': 'CelMek-small-body-extraction/1'})
                try:
                    with urllib.request.urlopen(request, timeout=180) as response:
                        raw = response.read()
                    payload = json.loads(raw)
                    if 'error' in payload or 'message' in payload:
                        raise ValueError('JPL error: ' + str(payload.get('error', payload.get('message'))))
                    signature = payload.get('signature', {})
                    expected = {SBDB_QUERY: '1.0', SBDB: '1.3', HORIZONS: '1.2'}[endpoint]
                    if not signature.get('version') or (expected and signature['version'] != expected):
                        raise ValueError('Unexpected API version; review schema: ' + str(signature))
                    envelope = dict(schema=1, url=url, fetched_at=timestamp(),
                                    response_sha256=hashlib.sha256(raw).hexdigest(), response=payload)
                    write_json(path, envelope)
                    break
                except urllib.error.HTTPError as error:
                    if error.code not in (429, 500, 502, 503, 504) or attempt == 3:
                        raise
                    time.sleep(2 ** attempt)
                except (urllib.error.URLError, TimeoutError):
                    if attempt == 3:
                        raise
                    time.sleep(2 ** attempt)
        signature = envelope['response'].get('signature', {})
        if signature.get('version') != {SBDB_QUERY: '1.0', SBDB: '1.3', HORIZONS: '1.2'}[endpoint]:
            raise ValueError('Unexpected cached API version; review schema: ' + str(signature))
        self.sources[key] = {k: envelope[k] for k in ('url', 'fetched_at', 'response_sha256')}
        return envelope['response'], key


def query_rows(client, parameters):
    payload, key = client.get(SBDB_QUERY, parameters)
    fields = payload.get('fields')
    data = payload.get('data', [])
    if not isinstance(fields, list) or not isinstance(data, list):
        raise ValueError('Invalid query response schema')
    rows = []
    for values in data:
        if len(values) != len(fields):
            raise ValueError('Query row width mismatch')
        row = dict(zip(fields, values))
        # Uncertainties can be asymmetric text, e.g. diameter_sigma='-1/+4'; preserve them.
        object_fields = {'spkid', 'pdes', 'full_name', 'name', 'prefix', 'class', 'pole', 'spec_B', 'spec_T', 'condition_code', 'orbit_id', 'diameter_sigma'}
        for field in fields:
            if field not in object_fields and row[field] is not None:
                number = float(row[field])
                if not math.isfinite(number):
                    raise ValueError('Nonfinite query field: ' + field)
                row[field] = number
        row['spkid'] = str(row['spkid'])
        if not row['spkid'].isdigit():
            raise ValueError('Invalid SPK identifier')
        row['source_cache_key'] = key
        rows.append(row)
    return rows


def numeric(row, field, default=math.inf):
    try:
        number = float(row[field])
        return number if math.isfinite(number) else default
    except (KeyError, TypeError, ValueError):
        return default


def select(client, args):
    manifest = json.loads(Path(args.manifest).read_text())
    if manifest.get('schema') != 1:
        raise ValueError('Unsupported selection manifest schema')
    bodies = {}
    def add(row, reason, kind):
        identifier = row['spkid']
        if identifier not in bodies:
            bodies[identifier] = dict(spkid=identifier, kind=kind, data=row, selection_reasons=[])
        if reason not in bodies[identifier]['selection_reasons']:
            bodies[identifier]['selection_reasons'].append(reason)

    asteroids = query_rows(client, {'fields': FIELDS, 'sb-kind': 'a', 'sort': '-diameter,spkid',
        'limit': str(args.asteroids), 'sb-cdata': json.dumps({'AND': ['diameter|DF']}, separators=(',', ':')),
        'full-prec': 'true'})
    if len(asteroids) != args.asteroids:
        raise ValueError('Insufficient asteroids with defined diameters')
    for rank, row in enumerate(asteroids, 1):
        add(row, f'Diameter ranking among cataloged values: {rank}', 'asteroid')
    for seed in manifest['mandatory_asteroids']:
        row = next((r for r in asteroids if r['pdes'] == seed['designation']), None)
        if row is None:
            # Fetch the exact designation even when no diameter is known.
            row_list = query_rows(client, {'fields': FIELDS, 'sb-kind': 'a', 'full-prec': 'true',
                'sb-cdata': json.dumps({'AND': ['pdes|EQ|' + seed['designation']]}, separators=(',', ':'))})
            if len(row_list) != 1:
                raise ValueError('Mandatory asteroid unresolved: ' + seed['designation'])
            row = row_list[0]
        add(row, seed['reason'], 'asteroid')
    if args.trojans:
        trojans = query_rows(client, {'fields': FIELDS, 'sb-kind': 'a', 'sb-class': 'TJN',
            'sort': '-diameter,spkid', 'limit': str(args.trojans), 'full-prec': 'true',
            'sb-cdata': json.dumps({'AND': ['diameter|DF']}, separators=(',', ':'))})
        if len(trojans) != args.trojans:
            raise ValueError('Insufficient Trojans with defined diameters')
        for rank, row in enumerate(trojans, 1):
            add(row, f'Jupiter Trojan diameter ranking: {rank}; L4/L5 camp unresolved', 'asteroid')

    candidates = query_rows(client, {'fields': COMET_FIELDS, 'sb-kind': 'c', 'sb-xfrag': 'true', 'full-prec': 'true'})
    selected = set()
    def add_comet(row, reason):
        if len(selected) < args.comets:
            selected.add(row['spkid']); add(row, reason, 'comet')
    missing_seeds = []
    for seed in manifest['comet_seeds']:
        if len(selected) >= args.comets:
            break
        rows = [r for r in candidates if r['pdes'] == seed['designation']]
        if len(rows) != 1:
            missing_seeds.append(seed['designation']); continue
        add_comet(rows[0], seed['reason'])
    if missing_seeds:
        raise ValueError('Curated comet seeds unresolved: ' + ', '.join(missing_seeds))
    # Fill remaining slots round-robin across orbit families. Within a family,
    # prefer known total magnitude M1, then known size, then primary designation.
    groups = {}
    for row in candidates:
        if row['spkid'] not in selected:
            groups.setdefault(row['class'], []).append(row)
    for group in groups.values():
        group.sort(key=lambda r: (numeric(r, 'M1'), -numeric(r, 'diameter', 0), r['pdes']))
    while len(selected) < args.comets:
        added = False
        for family in sorted(groups):
            if groups[family] and len(selected) < args.comets:
                row = groups[family].pop(0)
                add_comet(row, f'Representative orbit family {family}; ranked by catalog M1, then diameter')
                added = True
        if not added:
            raise ValueError('Insufficient comet candidates')
    result = dict(schema=1, generated_at=timestamp(), manifest=manifest,
        selection=dict(asteroids_by_known_diameter=args.asteroids, additional_trojan_quota=args.trojans,
                       comets=args.comets, comet_policy='Curated seeds then orbit-family diversity; not an objective interest ranking'),
        units=dict(diameter='km', GM='km^3/s^2', density='g/cm^3', rot_per='hours',
                   orbital_distances='AU', orbital_angles='degrees', epoch_and_tp='JD TDB'),
        sources=client.sources, bodies=sorted(bodies.values(), key=lambda b: (b['kind'], int(b['spkid']))))
    write_json(args.output, result)
    print(f"Wrote {len(result['bodies'])} POIs ({len(selected)} comets) to {args.output}")


def enrich(client, args):
    catalog = json.loads(Path(args.catalog).read_text())
    for body in catalog['bodies'][:args.limit]:
        payload, key = client.get(SBDB, {'spk': body['spkid'], 'phys-par': 'true', 'full-prec': 'true', 'alt-des': 'true'})
        if str(payload.get('object', {}).get('spkid')) != body['spkid']:
            raise ValueError('SBDB identity mismatch: ' + body['spkid'])
        body['sbdb_details'] = payload
        body['detail_source_cache_key'] = key
        print('Enriched', body['spkid'])
    catalog['sources'].update(client.sources)
    catalog['enriched_at'] = timestamp()
    catalog['enrichment'] = dict(enriched_count=sum('sbdb_details' in b for b in catalog['bodies']), total_count=len(catalog['bodies']))
    write_json(args.output, catalog)


def kernels(client, args):
    catalog = json.loads(Path(args.catalog).read_text())
    try:
        if dt.date.fromisoformat(args.start) >= dt.date.fromisoformat(args.stop):
            raise ValueError('start must precede stop')
    except ValueError as error:
        raise ValueError('Use ordered YYYY-MM-DD coverage dates') from error
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    existing = output/'manifest.json'
    previous = json.loads(existing.read_text()) if existing.exists() else {'schema': 1, 'kernels': [], 'sources': {}}
    if previous.get('schema') != 1 or any(k['requested_start_tdb'] != args.start or k['requested_stop_tdb'] != args.stop for k in previous['kernels']):
        raise ValueError('Kernel output directory has incompatible coverage; select a different directory')
    entries = {k['requested_spkid']: k for k in previous['kernels']}
    sources = previous['sources']
    for body in catalog['bodies'][:args.limit]:
        payload, key = client.get(HORIZONS, {'format': 'json', 'COMMAND': "'DES=" + body['spkid'] + ";'",
            'EPHEM_TYPE': 'SPK', 'OBJ_DATA': 'NO', 'START_TIME': "'" + args.start + "'",
            'STOP_TIME': "'" + args.stop + "'"})
        if 'spk' not in payload:
            raise ValueError('Horizons did not generate SPK for ' + body['spkid'] + ': ' + payload.get('result', ''))
        encoded = ''.join(payload['spk'].split())
        binary = base64.b64decode(encoded, validate=True)
        if not binary.startswith(b'DAF/SPK '):
            raise ValueError('Unexpected kernel header')
        identifier = str(payload.get('spk_file_id', ''))
        if not identifier.isdigit():
            raise ValueError('Invalid returned SPK identity')
        # Returned IDs can differ from requested aliases; preserve both.
        path = output / (identifier + '.bsp')
        temporary = path.with_suffix('.bsp.tmp'); temporary.write_bytes(binary); temporary.replace(path)
        entries[body['spkid']] = dict(requested_spkid=body['spkid'], returned_spkid=identifier, file=path.name,
            sha256=hashlib.sha256(binary).hexdigest(), requested_start_tdb=args.start,
            requested_stop_tdb=args.stop, source_cache_key=key)
        # Persist after each file so interrupted downloads are inspectable/resumable.
        sources.update(client.sources)
        write_json(output/'manifest.json', dict(schema=1, generated_at=timestamp(), sources=sources, kernels=list(entries.values())))
        print('Wrote', path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', default='Data/JPL/cache')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--refresh', action='store_true')
    commands = parser.add_subparsers(dest='command', required=True)
    selection = commands.add_parser('select')
    selection.add_argument('--manifest', default=str(DEFAULT_MANIFEST))
    selection.add_argument('--asteroids', type=int, default=1000)
    selection.add_argument('--trojans', type=int, default=50)
    selection.add_argument('--comets', type=int, default=100)
    selection.add_argument('--output', default='Data/JPL/pois.json')
    for name in ('enrich', 'kernels'):
        command = commands.add_parser(name)
        command.add_argument('--catalog', default='Data/JPL/pois.json')
        command.add_argument('--limit', type=int)
        command.add_argument('--output', default='Data/JPL/pois-enriched.json' if name == 'enrich' else 'Data/JPL/SPK')
        if name == 'kernels':
            command.add_argument('--start', required=True)
            command.add_argument('--stop', required=True)
    args = parser.parse_args()
    if args.offline and args.refresh:
        parser.error('--offline and --refresh conflict')
    if args.command == 'select' and (args.asteroids < 1 or args.comets < 1 or args.trojans < 0):
        parser.error('positive asteroid/comet counts and nonnegative Trojan quota required')
    if getattr(args, 'limit', None) is not None and args.limit < 1:
        parser.error('--limit must be positive')
    client = Client(args.cache, args.offline, args.refresh)
    globals()[args.command](client, args)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        print('Extraction failed:', error, file=sys.stderr)
        sys.exit(1)
