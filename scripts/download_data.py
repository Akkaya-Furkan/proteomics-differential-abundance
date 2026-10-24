"""Download only public PRIDE data; verify pinned bytes before atomic replacement."""
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://ftp.ebi.ac.uk/pride-archive/2014/09/PXD000279/proteomebenchmark.zip'
GO_URL = 'https://rest.uniprot.org/uniprotkb/stream?query=proteome%3AUP000000625&format=tsv&fields=accession,go_p'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def download(url, path, expected):
    path = Path(path)
    if path.exists() and sha256(path) == expected:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.part')
    try:
        with urllib.request.urlopen(url, timeout=180) as response, tmp.open('wb') as out:
            while chunk := response.read(1024 * 1024):
                out.write(chunk)
        if sha256(tmp) != expected:
            raise ValueError(f'Checksum mismatch: {url}; destination unchanged')
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def main():
    manifest = json.loads((ROOT / 'data/sources.json').read_text())
    for item in manifest['sources']:
        download(item['url'], ROOT / item['path'], item['sha256'])
    archive = ROOT / 'data/raw/proteomebenchmark.zip'
    with zipfile.ZipFile(archive) as z:
        for name, digest in manifest['extracted'].items():
            content = z.read(name)  # Explicit allowlist; no extractall/path traversal.
            if hashlib.sha256(content).hexdigest() != digest:
                raise ValueError(f'Archive member checksum mismatch: {name}')
        for name in manifest['extracted']:
            target = ROOT / 'data/raw' / name
            tmp = target.with_suffix('.part')
            tmp.write_bytes(z.read(name))
            tmp.replace(target)
    print('Verified 2/2 public sources and 2/2 extracted files.')


if __name__ == '__main__':
    main()
