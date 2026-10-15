"""Serve a built DocC archive locally so its JavaScript and JSON assets can load.
Usage: python3 Tools/Documentation/preview.py [--port 8080] [--archive PATH]
"""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit


class DocumentationHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        route = urlsplit(self.path).path
        if route == "/":
            self.send_response(302)
            self.send_header("Location", "/documentation")
            self.end_headers()
            return
        # DocC uses client-side routes. Serve the corresponding static entry point
        # without SimpleHTTPRequestHandler's automatic trailing-slash redirect.
        if route == "/documentation" or route.startswith("/documentation/") or route.startswith("/tutorials/"):
            path = Path(self.translate_path(route))
            if path.is_dir() and (path / "index.html").is_file():
                self.path = route.rstrip("/") + "/index.html"
            elif not path.exists() and not path.suffix:
                self.path = "/index.html"
        super().do_GET()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    parser.add_argument("--archive", type=Path, default=root / ".build/documentation/Physics.doccarchive")
    parser.add_argument("--port", type=int, default=8080)
    args = parser.parse_args()
    archive = args.archive.resolve()
    if not (archive / "index.html").is_file() or not (archive / "data/documentation.json").is_file():
        parser.error(f"No combined DocC archive at {archive}. Run Tools/Documentation/build.py first.")
    handler = partial(DocumentationHandler, directory=str(archive))
    with ThreadingHTTPServer(("127.0.0.1", args.port), handler) as server:
        print(f"Documentation: http://127.0.0.1:{server.server_port}/documentation", flush=True)
        print("Press Ctrl-C to stop the preview server.", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
