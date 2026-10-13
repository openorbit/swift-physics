# SPDX-License-Identifier: Apache-2.0
import base64
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch
import urllib.parse
import extract

class FakeClient:
    def __init__(self, payload):
        self.payload = payload
        self.sources = {}
    def get(self, endpoint, parameters):
        self.sources['fake'] = {'url': endpoint, 'fetched_at': 'fixture', 'response_sha256': 'fixture'}
        return self.payload, 'fake'

class ExtractionTests(unittest.TestCase):
    def test_query_normalizes_numbers_preserves_asymmetric_uncertainty_and_unknowns(self):
        client = FakeClient({'fields': ['spkid','diameter','diameter_sigma','GM'],
                             'data': [[20000001, '10.2', '-1/+4', None]]})
        row = extract.query_rows(client, {})[0]
        self.assertEqual(row['spkid'], '20000001')
        self.assertEqual(row['diameter'], 10.2)
        self.assertEqual(row['diameter_sigma'], '-1/+4')
        self.assertIsNone(row['GM'])

    def test_query_rejects_nonfinite_numbers_and_wrong_row_width(self):
        for payload in [{'fields':['spkid','diameter'],'data':[[1,'NaN']]},
                        {'fields':['spkid','diameter'],'data':[[1]]}]:
            with self.assertRaises(ValueError):
                extract.query_rows(FakeClient(payload), {})

    def test_offline_cache_hit_miss_and_schema_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            parameters = {'fields': 'spkid'}
            url = extract.SBDB_QUERY + '?' + urllib.parse.urlencode(sorted(parameters.items()))
            key = hashlib.sha256(url.encode()).hexdigest()
            path = Path(directory)/(key+'.json')
            envelope = dict(schema=1,url=url,fetched_at='fixture',response_sha256='hash',
                            response={'signature':{'version':'1.0'},'fields':['spkid'],'data':[[1]]})
            extract.write_json(path,envelope)
            client = extract.Client(directory, offline=True)
            self.assertEqual(client.get(extract.SBDB_QUERY,parameters)[1], key)
            with self.assertRaisesRegex(ValueError,'cache miss'):
                client.get(extract.SBDB_QUERY,{'fields':'name'})
            envelope['response']['signature']['version']='99'
            extract.write_json(path,envelope)
            with self.assertRaisesRegex(ValueError,'version'):
                client.get(extract.SBDB_QUERY,parameters)

    def test_selection_deduplicates_mandatory_and_trojan_and_fills_families(self):
        def row(id, pdes, family='MBA', magnitude=None):
            return dict(spkid=str(id),pdes=pdes,**{'class':family},M1=magnitude,diameter=10)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            manifest=dict(schema=1,mandatory_asteroids=[dict(designation='1',reason='mandatory')],
                          comet_seeds=[dict(designation='1P',reason='curated')])
            extract.write_json(root/'selection.json',manifest)
            args=SimpleNamespace(manifest=root/'selection.json',asteroids=2,trojans=1,comets=3,output=root/'pois.json')
            responses=[[row(1,'1'),row(2,'2')],[row(2,'2','TJN')],
                       [row(11,'1P','JFC'),row(12,'2P','JFC',3),row(13,'x','HYP',1),row(14,'y','JFC',None)]]
            with patch.object(extract,'query_rows',side_effect=responses):
                extract.select(FakeClient({}),args)
            catalog=json.loads(args.output.read_text())
            self.assertEqual(len(catalog['bodies']),5)
            self.assertEqual(len(catalog['bodies'][0]['selection_reasons']),2)
            self.assertEqual(len(catalog['bodies'][1]['selection_reasons']),2)
            self.assertEqual({b['spkid'] for b in catalog['bodies'] if b['kind']=='comet'}, {'11','12','13'})

    def test_kernel_decode_manifest_resume_and_coverage_conflicts(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            extract.write_json(root/'catalog.json',{'bodies':[{'spkid':'1'},{'spkid':'2'}]})
            args=SimpleNamespace(catalog=root/'catalog.json',output=root/'kernels',limit=1,start='2000-01-01',stop='2000-01-03')
            payload={'spk':base64.b64encode(b'DAF/SPK '+b'x'*1016).decode(),'spk_file_id':'1'}
            extract.kernels(FakeClient(payload),args)
            self.assertEqual((args.output/'1.bsp').read_bytes()[:8],b'DAF/SPK ')
            args.limit=2
            class Sequential(FakeClient):
                def get(self,endpoint,parameters):
                    self.payload['spk_file_id']='1' if 'DES=1;' in parameters['COMMAND'] else '2'
                    return super().get(endpoint,parameters)
            extract.kernels(Sequential(payload),args)
            self.assertEqual(len(json.loads((args.output/'manifest.json').read_text())['kernels']),2)
            args.stop='2001-01-01'
            with self.assertRaisesRegex(ValueError,'incompatible coverage'):
                extract.kernels(FakeClient(payload),args)

    def test_kernel_rejects_missing_spk_bad_base64_and_wrong_header(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            extract.write_json(root/'catalog.json',{'bodies':[{'spkid':'1'}]})
            args=SimpleNamespace(catalog=root/'catalog.json',output=root/'kernels',limit=None,start='2000-01-01',stop='2000-01-03')
            for payload in [{}, {'spk':'!!!'}, {'spk':base64.b64encode(b'bad').decode()}]:
                with self.assertRaises(ValueError):
                    extract.kernels(FakeClient(payload),args)
            self.assertFalse((args.output/'manifest.json').exists())

    def test_enrichment_rejects_identity_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            extract.write_json(root/'catalog.json',{'bodies':[{'spkid':'1'}]})
            args=SimpleNamespace(catalog=root/'catalog.json',output=root/'output.json',limit=None)
            with self.assertRaisesRegex(ValueError,'identity mismatch'):
                extract.enrich(FakeClient({'object':{'spkid':'2'}}),args)
            self.assertFalse(args.output.exists())

if __name__ == '__main__':
    unittest.main()
