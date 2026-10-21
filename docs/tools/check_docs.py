#!/usr/bin/env python3
"""Check local HTML links, anchors, assets and authored code block wrappers."""
import argparse
import ast
import hashlib
import json
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

class Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids, self.refs, self.pre_depth, self.pre_wrapped = set(), [], 0, []
        self.current_count = 0
    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs: self.ids.add(attrs['id'])
        if tag == 'a' and 'name' in attrs: self.ids.add(attrs['name'])
        if tag == 'a' and 'aria-current' in attrs and attrs['aria-current'] == 'page': self.current_count += 1
        if tag in ('a', 'link', 'script', 'img'):
            key = 'href' if tag in ('a', 'link') else 'src'
            if attrs.get(key): self.refs.append((tag, attrs[key]))
        if tag == 'pre':
            self.pre_depth += 1
            self.pre_wrapped.append('code-block' in attrs.get('class', '').split())
    def handle_endtag(self, tag):
        if tag == 'pre' and self.pre_depth:
            self.pre_depth -= 1

def inspect(html):
    page = Page(); page.feed(html.read_text(errors='replace')); return page

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('site', type=Path, help='built docs/site directory')
    args = parser.parse_args()
    site = args.site.resolve()
    if not (site / 'index.html').is_file() or not (site / 'api/index.html').is_file():
        parser.error('expected built site/index.html and site/api/index.html')
    html_files = sorted(site.rglob('*.html'))
    parsed = {p: inspect(p) for p in html_files}
    errors = []
    required_pages = {
        'index.html', 'installation.html', 'first-seismogram.html',
        'python.html', 'plotting.html', 'kernels.html',
        'library-example.html', 'examples.html', 'input-files.html',
        'model-and-output-formats.html', 'methods-and-validation.html',
        'release-notes.html', 'maintainers.html',
    }
    missing_pages = sorted(required_pages - {p.relative_to(site).as_posix() for p in html_files})
    if missing_pages:
        errors.append('missing authored pages: ' + ', '.join(missing_pages))
    required_api = {
        'api/Input_8hpp.html', 'api/Instrument_8hpp.html',
        'api/Run_8hpp.html', 'api/Synthetics_8hpp.html',
        'api/Kernels_8hpp.html', 'api/classYSpec_1_1Model.html',
        'api/classYSpec_1_1RadialSolver.html',
    }
    authored_refs = {
        ref for path, page in parsed.items()
        if 'api' not in path.relative_to(site).parts
        for tag, ref in page.refs if tag == 'a'
    }
    missing_api = sorted(required_api - {urlsplit(ref).path for ref in authored_refs})
    if missing_api:
        errors.append('authored guides must link generated API pages: ' + ', '.join(missing_api))
    required_api_fragments = {
        'api/Input_8hpp.html#a159b2668249241093472863aae07800b',
        'api/Instrument_8hpp.html#ae4d3fd526efb1cdb3c82d9965a3f22d4',
        'api/Run_8hpp.html#a005fd5bf5085a7ce92b77071a426df19',
        'api/Synthetics_8hpp.html#a37f382b0fc276cb40ecc2f2a948e101c',
        'api/Kernels_8hpp.html#a87ca499eeb4cdae1c1060626cd687760',
        'api/Input_8hpp.html#a9f0b0ad165807071cecf1808d796bca3',
        'api/Run_8hpp.html#a94dff6e2a13c21656f02548b0db61ba2',
        'api/Kernels_8hpp.html#a2700b3e83da66c395e2b58f31a290f40',
        'api/classYSpec_1_1Model.html#a11f4503d9302dc5a1f79020f2cc6b72b',
    }
    actual_api_fragments = {
        f'{urlsplit(ref).path}#{urlsplit(ref).fragment}'
        for ref in authored_refs if urlsplit(ref).fragment
    }
    missing_fragments = sorted(required_api_fragments - actual_api_fragments)
    if missing_fragments:
        errors.append('authored guides must link exact generated API symbols: ' + ', '.join(missing_fragments))
    for page_path, page in parsed.items():
        for tag, ref in page.refs:
            url = urlsplit(ref)
            if url.scheme or url.netloc or ref.startswith('//'):
                continue
            target = (page_path.parent / unquote(url.path)).resolve() if url.path else page_path
            try: target.relative_to(site)
            except ValueError:
                errors.append(f'{page_path.relative_to(site)}: local URL escapes site: {ref}')
                continue
            if not target.is_file():
                errors.append(f'{page_path.relative_to(site)}: missing {tag} target {ref}')
                continue
            fragment = unquote(url.fragment)
            if fragment and target.suffix.lower() == '.html':
                target_page = parsed.get(target) or inspect(target)
                if fragment not in target_page.ids:
                    errors.append(f'{page_path.relative_to(site)}: missing anchor {ref}')
        if 'api' not in page_path.relative_to(site).parts:
            rel = page_path.relative_to(site)
            source = page_path.read_text(errors='replace')
            if '{{' in source or '}}' in source:
                errors.append(f'{rel}: unresolved build placeholder remains')
            if page.current_count != 1:
                errors.append(f'{rel}: expected exactly one aria-current="page" navigation link')
            if page.pre_wrapped and not all(page.pre_wrapped):
                errors.append(f'{rel}: each authored pre block must use the copy-code wrapper')
            if not (site / 'assets/site.js').is_file():
                errors.append(f'{rel}: shared code-copy script is missing')
    js = (site / 'assets/site.js').read_text()
    for behavior in ('navigator.clipboard.writeText', 'document.execCommand', 'aria-live'):
        if behavior not in js: errors.append(f'copy script lacks {behavior}')
    try:
        json.loads((site / 'examples/first-seismogram.json').read_text())
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f'tutorial JSON fixture is invalid or missing: {error}')
    if not (site / 'assets/first-seismogram.svg').is_file():
        errors.append('tutorial plot asset is missing')
    try:
        json.loads((site / 'examples/snapshot-small.json').read_text())
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f'snapshot JSON fixture is invalid or missing: {error}')
    try:
        json.loads((site / 'examples/workflows.json').read_text())
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f'workflow JSON fixture is invalid or missing: {error}')
    for fixture in ('kernels-small.json', 'kernels-volume-small.json'):
        try:
            json.loads((site / 'examples' / fixture).read_text())
        except (OSError, json.JSONDecodeError) as error:
            errors.append(f'{fixture} is invalid or missing: {error}')
    python_example = site / 'examples/python-small-calculation.py'
    try:
        ast.parse(python_example.read_text(), filename=str(python_example))
    except (OSError, SyntaxError) as error:
        errors.append(f'Python calculation example is invalid or missing: {error}')
    for required_asset in ('examples/prem.200', 'examples/station.pz',
                           'examples/workflows.json',
                           'examples/library-workflows.cpp',
                           'examples/kernels-small.json',
                           'examples/kernels-volume-small.json',
                           'examples/python-small-calculation.py',
                           'assets/traces-small.png',
                           'assets/kernels-radial-small.png',
                           'assets/kernels-section-small.png',
                           'assets/adjoint-kernels.pdf'):
        if not (site / required_asset).is_file():
            errors.append(f'required downloadable asset is missing: {required_asset}')
    try:
        fixture_hash = hashlib.sha256((site / 'examples/prem.200').read_bytes()).digest()
        canonical_hash = hashlib.sha256((Path(__file__).resolve().parents[2] / 'examples/prem.200').read_bytes()).digest()
        if fixture_hash != canonical_hash:
            errors.append('published PREM model copy differs from examples/prem.200')
    except OSError as error:
        errors.append(f'could not verify the published PREM model copy: {error}')
    if errors:
        print('\n'.join(errors), file=sys.stderr); return 1
    print(f'Checked {len(html_files)} HTML pages, local links, anchors, assets and code-copy wrappers.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
