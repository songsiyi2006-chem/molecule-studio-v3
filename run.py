"""Start the local app; optionally open the browser after health becomes ready."""
import argparse
from pathlib import Path
import os
import socket
import threading
import time
import urllib.request
import webbrowser


def open_when_ready(url):
    for _ in range(120):
        try:
            with urllib.request.urlopen(url + "/health", timeout=1) as r:
                if r.status == 200:
                    webbrowser.open(url)
                    return
        except OSError:
            pass
        time.sleep(0.5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Molecule Studio V3 local launcher")
    parser.add_argument("--port", type=int, default=8003)
    parser.add_argument("--open", action="store_true", help="Open browser when ready")
    args = parser.parse_args()
    os.chdir(Path(__file__).resolve().parent)
    if not 1 <= args.port <= 65535:
        raise SystemExit("Port must be between 1 and 65535.")
    try:
        with socket.socket() as s:
            s.bind(("127.0.0.1", args.port))
    except OSError:
        raise SystemExit(f"Port {args.port} is in use. Try: python run.py --port {args.port + 1} --open")
    for variable in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        os.environ.setdefault(variable, '2')
    import uvicorn
    if args.open:
        threading.Thread(target=open_when_ready, args=(f"http://127.0.0.1:{args.port}",), daemon=True).start()
    uvicorn.run("app.main:app", host="127.0.0.1", port=args.port, log_level="info")
