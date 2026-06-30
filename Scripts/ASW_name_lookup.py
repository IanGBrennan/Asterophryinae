#!/usr/bin/env python3
"""
ASW Species Name Updater
========================
Looks up old amphibian species names against the Amphibian Species of the
World (ASW) database (https://amphibiansoftheworld.amnh.org/) using the
site's BASIC SEARCH feature, and returns the current accepted name.

HOW IT WORKS
------------
For each input name (e.g. "Barygenys maculata") the script queries ASW's
basic search using an exact quoted phrase:

    https://amphibiansoftheworld.amnh.org/amphib/basic_search?basic_query="Barygenys+maculata"

This mirrors what you would type into the Basic Search box on the ASW
website.  Because the phrase is quoted, it only returns records that
contain that exact string — typically the species page where the name
is the accepted name, or (if it's a synonym) the page where it appears
in the synonymy section.

The first result link pointing to an ASW species page is taken as the
answer.  The accepted name is read from the page title.

FALLBACK
--------
If the basic search returns no hits, the script retries without quotes
so a looser match is attempted.  These rows are flagged with status
'fallback_unquoted' in the output so you know they need manual checking.

Usage:
    pip install requests beautifulsoup4
    python asw_name_lookup.py input.csv output.csv

Options:
    --column COLUMN   Column containing species names (auto-detected if omitted)
    --delay  SECONDS  Pause between requests in seconds (default 1.5)

Example:
    python asw_name_lookup.py frogs.csv frogs_updated.csv
    python asw_name_lookup.py frogs.csv frogs_updated.csv --column old_name --delay 2
"""

import csv
import sys
import time
import argparse
import re
import urllib.parse
from pathlib import Path

try:
    import requests
    from bs4 import BeautifulSoup
except ImportError:
    print("ERROR: Missing dependencies. Please run:")
    print("  pip install requests beautifulsoup4")
    sys.exit(1)


# ── Constants ──────────────────────────────────────────────────────────────────

ASW_BASE      = "https://amphibiansoftheworld.amnh.org"
ASW_SEARCH    = ASW_BASE + "/amphib/basic_search"

# Mimic a regular browser so the server doesn't block us
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0.0.0 Safari/537.36"
    ),
    "Accept":          "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
    "Referer":         ASW_BASE,
}

# Column names to auto-detect for species names (case-insensitive)
AUTO_DETECT_COLS = {"species", "species_name", "old_name", "name", "taxon", "old_species"}

# Regex matching ASW species-page URLs:
#   /Amphibia/<order>/<family>/.../<Genus>/<Genus>-<species>
SPECIES_URL_RE = re.compile(r"/Amphibia/.+/[A-Z][a-z]+-[a-z]")


# ── HTTP helper ────────────────────────────────────────────────────────────────

def get(session: requests.Session, url: str) -> requests.Response | None:
    """GET with error handling; returns None on any non-200 or network error."""
    try:
        r = session.get(url, headers=HEADERS, timeout=20, allow_redirects=True)
        return r if r.status_code == 200 else None
    except requests.RequestException:
        return None


# ── Page parsing ───────────────────────────────────────────────────────────────

def is_species_page(soup: BeautifulSoup) -> bool:
    """True if the page looks like an ASW species record (not a search-results page)."""
    title = soup.find("title")
    if not title:
        return False
    # Species page titles: "Genus species (Author, Year) | Amphibian Species of the World"
    return bool(re.match(r"^[A-Z][a-z]+ [a-z]+", title.get_text(strip=True)))


def extract_accepted_name(soup: BeautifulSoup, page_url: str) -> str:
    """
    Pull the accepted binomial from an ASW species page.
    Tries the <title> tag first, then the page URL slug as a fallback.
    """
    title = soup.find("title")
    if title:
        m = re.match(r"^([A-Z][a-z]+(?: [a-z]+){1,2})", title.get_text(strip=True))
        if m:
            return m.group(1).strip()

    # Fallback: decode from URL — last two path segments are /Genus/Genus-species
    m = re.search(r"/([A-Z][a-z]+)/([A-Z][a-z]+-[a-z]+(?:-[a-z]+)?)/?$", page_url)
    if m:
        return m.group(2).replace("-", " ")

    return ""


def first_species_link(soup: BeautifulSoup) -> str | None:
    """
    From a basic-search results page, return the href of the first link
    that points to an ASW species record page.
    """
    for a in soup.find_all("a", href=True):
        if SPECIES_URL_RE.search(a["href"]):
            href = a["href"]
            return href if href.startswith("http") else ASW_BASE + href
    return None


# ── Core lookup ────────────────────────────────────────────────────────────────

def basic_search(query: str, session: requests.Session, quoted: bool = True) -> requests.Response | None:
    """
    Hit the ASW basic search endpoint.
    When quoted=True  the query is wrapped in double-quotes for an exact phrase match.
    When quoted=False a plain (unquoted) keyword search is used as a fallback.
    """
    phrase = f'"{query}"' if quoted else query
    params = urllib.parse.urlencode({"basic_query": phrase})
    url = f"{ASW_SEARCH}?{params}"
    return get(session, url)


def search_asw(species_name: str, session: requests.Session) -> dict:
    """
    Look up one species name via ASW basic search.

    Returns a dict with:
        accepted_name  – current ASW accepted binomial
        accepted_url   – URL of the ASW species page
        status         – 'accepted' | 'synonym' | 'fallback_unquoted' |
                         'not_found' | 'error'
        notes          – human-readable detail
    """
    empty = {"accepted_name": "", "accepted_url": "", "status": "", "notes": ""}

    # ── Step 1: exact quoted phrase search ────────────────────────────────────
    r = basic_search(species_name, session, quoted=True)

    if r is None:
        return {**empty, "status": "error",
                "notes": "Basic search request failed (network error or server blocked)"}

    soup = BeautifulSoup(r.text, "html.parser")

    # The search may redirect straight to a species page (single exact hit)
    if is_species_page(soup):
        name = extract_accepted_name(soup, r.url)
        return _build_result(name, r.url, species_name)

    # Otherwise parse the search-results listing
    sp_url = first_species_link(soup)

    if sp_url:
        sp_r = get(session, sp_url)
        if sp_r:
            sp_soup = BeautifulSoup(sp_r.text, "html.parser")
            name = extract_accepted_name(sp_soup, sp_r.url)
            return _build_result(name, sp_r.url, species_name)

    # ── Step 2: fallback — unquoted search ───────────────────────────────────
    r2 = basic_search(species_name, session, quoted=False)

    if r2 is None:
        return {**empty, "status": "not_found",
                "notes": "No results found in quoted or unquoted basic search"}

    soup2 = BeautifulSoup(r2.text, "html.parser")

    if is_species_page(soup2):
        name = extract_accepted_name(soup2, r2.url)
        return {
            "accepted_name": name,
            "accepted_url":  r2.url,
            "status":        "fallback_unquoted",
            "notes": (
                "⚠ Not found by exact quoted search. Result is from an unquoted fallback "
                "search and may not be correct — please verify manually."
            ),
        }

    sp_url2 = first_species_link(soup2)
    if sp_url2:
        sp_r2 = get(session, sp_url2)
        if sp_r2:
            sp_soup2 = BeautifulSoup(sp_r2.text, "html.parser")
            name = extract_accepted_name(sp_soup2, sp_r2.url)
            return {
                "accepted_name": name,
                "accepted_url":  sp_r2.url,
                "status":        "fallback_unquoted",
                "notes": (
                    "⚠ Not found by exact quoted search. Result is from an unquoted fallback "
                    "search and may not be correct — please verify manually at: " + sp_r2.url
                ),
            }

    return {**empty, "status": "not_found",
            "notes": "No results found in ASW basic search (quoted or unquoted)"}


def _build_result(accepted_name: str, page_url: str, input_name: str) -> dict:
    """Classify the result and return the standard result dict."""
    if not accepted_name:
        return {
            "accepted_name": "",
            "accepted_url":  page_url,
            "status":        "error",
            "notes":         "Species page found but accepted name could not be parsed",
        }

    input_parts   = input_name.strip().lower().split()
    accepted_parts = accepted_name.strip().lower().split()

    # If at least genus + epithet match, the name is already the accepted name
    if (len(input_parts) >= 2 and len(accepted_parts) >= 2
            and input_parts[0] == accepted_parts[0]
            and input_parts[1] == accepted_parts[1]):
        status = "accepted"
        notes  = "Name matches current accepted name in ASW"
    else:
        status = "synonym"
        notes  = (
            f"Input name '{input_name}' found in synonymy; "
            f"current accepted name is '{accepted_name}'"
        )

    return {
        "accepted_name": accepted_name,
        "accepted_url":  page_url,
        "status":        status,
        "notes":         notes,
    }


# ── CSV handling ───────────────────────────────────────────────────────────────

def detect_name_column(fieldnames: list[str], override: str | None) -> str:
    if override:
        if override not in fieldnames:
            print(f"ERROR: Column '{override}' not found in CSV.")
            print(f"  Available columns: {', '.join(fieldnames)}")
            sys.exit(1)
        return override
    for col in fieldnames:
        if col.strip().lower() in AUTO_DETECT_COLS:
            return col
    print(f"WARNING: Could not auto-detect species column. Using first column: '{fieldnames[0]}'")
    print(f"  Use --column to specify the correct column.")
    return fieldnames[0]


def process_csv(input_path: Path, output_path: Path, column: str | None, delay: float):
    session = requests.Session()

    with open(input_path, newline="", encoding="utf-8-sig") as f_in:
        reader = csv.DictReader(f_in)
        if not reader.fieldnames:
            print("ERROR: Input CSV appears to be empty.")
            sys.exit(1)
        name_col  = detect_name_column(list(reader.fieldnames), column)
        new_fields = ["asw_accepted_name", "asw_status", "asw_url", "asw_notes"]
        out_fields = list(reader.fieldnames) + new_fields
        rows       = list(reader)

    print(f"Using column : '{name_col}'")
    print(f"Input rows   : {len(rows)}")
    print()

    STATUS_LABEL = {
        "accepted":           "ACCEPTED",
        "synonym":            "SYNONYM →",
        "fallback_unquoted":  "FALLBACK ⚠",
        "not_found":          "NOT FOUND",
        "error":              "ERROR",
        "skipped":            "SKIPPED",
    }

    results = []
    for i, row in enumerate(rows, 1):
        species = row.get(name_col, "").strip()
        if not species:
            print(f"  [{i:>4}/{len(rows)}] SKIPPED (empty name)")
            row.update({"asw_accepted_name": "", "asw_status": "skipped",
                        "asw_url": "", "asw_notes": "Empty name"})
            results.append(row)
            continue

        print(f"  [{i:>4}/{len(rows)}] {species}", end=" ... ", flush=True)
        lookup = search_asw(species, session)

        label   = STATUS_LABEL.get(lookup["status"], lookup["status"].upper())
        display = lookup["accepted_name"] or "(not found)"
        suffix  = f" {display}" if lookup["status"] == "synonym" else f"  {display}"
        print(f"{label}{suffix}")

        row["asw_accepted_name"] = lookup["accepted_name"]
        row["asw_status"]        = lookup["status"]
        row["asw_url"]           = lookup["accepted_url"]
        row["asw_notes"]         = lookup["notes"]
        results.append(row)

        if i < len(rows):
            time.sleep(delay)

    with open(output_path, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(results)

    statuses = [r["asw_status"] for r in results]
    n_fb     = statuses.count("fallback_unquoted")

    print(f"\n{'='*60}")
    print(f"Output written to : {output_path}")
    print(f"  Accepted (name unchanged)  : {statuses.count('accepted')}")
    print(f"  Synonym  (name updated)    : {statuses.count('synonym')}")
    print(f"  Fallback (verify manually) : {n_fb}")
    print(f"  Not found                  : {statuses.count('not_found')}")
    print(f"  Error                      : {statuses.count('error')}")
    print(f"  Skipped  (empty)           : {statuses.count('skipped')}")
    if n_fb:
        print()
        print("  ⚠  Fallback rows were found via an unquoted search only.")
        print("     They may match a different taxon. Check asw_url to verify.")
    print(f"{'='*60}")
    print()
    print("Citation:")
    print("  Frost, D.R. 2026. Amphibian Species of the World: an Online Reference.")
    print("  Version 6.2. AMNH. https://amphibiansoftheworld.amnh.org/index.php")


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Update amphibian species names via the ASW basic search.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input",   type=Path, help="Input CSV with old species names")
    parser.add_argument("output",  type=Path, help="Output CSV with updated names")
    parser.add_argument("--column", default=None,
                        help="Column containing species names (auto-detected if omitted)")
    parser.add_argument("--delay",  type=float, default=1.5,
                        help="Seconds to wait between requests (default: 1.5)")
    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)

    process_csv(args.input, args.output, args.column, args.delay)


if __name__ == "__main__":
    main()
