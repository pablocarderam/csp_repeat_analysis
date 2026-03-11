#!/usr/bin/env python3
"""
Parse GenBank peptide sequence file and extract accession, sequence, and geo_loc_name to CSV.

This code was generated using Anthropic's Claude.
"""

import csv
import sys

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError


def classify_region(geo_loc_name, isolate_name="", reference_title="", journal=""):
    """
    Classify geographic location into regions.

    Args:
        geo_loc_name: Geographic location string
        isolate_name: Isolate tag string (fallback for classification)
        reference_title: First reference title (second fallback for classification)
        journal: Journal field from first reference (third fallback for classification)

    Returns:
        Tuple of (region classification, source used: 'geo', 'isolate', 'reference', 'journal', or '')
    """
    if not geo_loc_name and not isolate_name and not reference_title and not journal:
        return "", ""

    # Try geo_loc_name first
    location_lower = geo_loc_name.lower() if geo_loc_name else ""

    # Africa
    africa_countries = [
        "algeria",
        "angola",
        "benin",
        "botswana",
        "burkina faso",
        "burundi",
        "cameroon",
        "cape verde",
        "central african republic",
        "chad",
        "comoros",
        "congo",
        "djibouti",
        "egypt",
        "equatorial guinea",
        "eritrea",
        "ethiopia",
        "gabon",
        "gambia",
        "ghana",
        "guinea",
        "guinea-bissau",
        "ivory coast",
        "kenya",
        "lesotho",
        "liberia",
        "libya",
        "madagascar",
        "malawi",
        "mali",
        "mauritania",
        "mauritius",
        "morocco",
        "mozambique",
        "namibia",
        "niger",
        "nigeria",
        "rwanda",
        "senegal",
        "seychelles",
        "sierra leone",
        "somalia",
        "south africa",
        "south sudan",
        "sudan",
        "tanzania",
        "togo",
        "tunisia",
        "uganda",
        "zambia",
        "zimbabwe",
    ]

    # Add major African cities and regions
    africa_locations = africa_countries + [
        "yaoundé",
        "yaounde",
        "nairobi",
        "lagos",
        "cairo",
        "kinshasa",
        "luanda",
        "accra",
        "khartoum",
        "dakar",
        "abidjan",
        "cape town",
        "johannesburg",
        "casablanca",
        "addis ababa",
        "kampala",
        "dar es salaam",
        "asembo",
    ]

    # West Asia (Middle East)
    west_asia_countries = [
        "afghanistan",
        "armenia",
        "azerbaijan",
        "bahrain",
        "cyprus",
        "georgia",
        "iran",
        "iraq",
        "israel",
        "jordan",
        "kuwait",
        "lebanon",
        "oman",
        "palestine",
        "qatar",
        "saudi arabia",
        "syria",
        "turkey",
        "united arab emirates",
        "uae",
        "yemen",
    ]

    west_asia_locations = west_asia_countries + [
        "tehran",
        "baghdad",
        "damascus",
        "beirut",
        "amman",
        "riyadh",
        "dubai",
        "abu dhabi",
        "doha",
        "kuwait city",
        "ankara",
        "istanbul",
        "baluchistan",
        "balochistan",
        "baluchistan 1",
    ]

    # South Asia
    south_asia_countries = [
        "bangladesh",
        "bhutan",
        "india",
        "maldives",
        "nepal",
        "pakistan",
        "sri lanka",
    ]

    south_asia_locations = south_asia_countries + [
        "delhi",
        "mumbai",
        "kolkata",
        "chennai",
        "bangalore",
        "hyderabad",
        "dhaka",
        "karachi",
        "lahore",
        "islamabad",
        "colombo",
        "kathmandu",
    ]

    # Southeast Asia
    southeast_asia_countries = [
        "brunei",
        "cambodia",
        "indonesia",
        "laos",
        "malaysia",
        "myanmar",
        "burma",
        "philippines",
        "singapore",
        "thailand",
        "timor-leste",
        "east timor",
        "vietnam",
        "viet nam",
    ]

    southeast_asia_locations = southeast_asia_countries + [
        "bangkok",
        "jakarta",
        "manila",
        "hanoi",
        "ho chi minh",
        "kuala lumpur",
        "yangon",
        "phnom penh",
        "vientiane",
        "thai",
    ]

    # East Asia
    east_asia_countries = [
        "china",
        "hong kong",
        "japan",
        "mongolia",
        "north korea",
        "south korea",
        "korea",
        "taiwan",
        "macau",
    ]

    east_asia_locations = east_asia_countries + [
        "beijing",
        "shanghai",
        "guangzhou",
        "shenzhen",
        "tokyo",
        "osaka",
        "seoul",
        "taipei",
        "ulaanbaatar",
    ]

    # Americas (North, Central, and South America)
    americas_countries = [
        # North America
        "canada",
        "united states",
        "usa",
        "mexico",
        # Central America and Caribbean
        "belize",
        "costa rica",
        "el salvador",
        "guatemala",
        "honduras",
        "nicaragua",
        "panama",
        "antigua",
        "barbuda",
        "bahamas",
        "barbados",
        "cuba",
        "dominica",
        "dominican republic",
        "grenada",
        "haiti",
        "jamaica",
        "saint kitts",
        "nevis",
        "saint lucia",
        "saint vincent",
        "grenadines",
        "trinidad",
        "tobago",
        # South America
        "argentina",
        "bolivia",
        "brazil",
        "chile",
        "colombia",
        "ecuador",
        "guyana",
        "paraguay",
        "peru",
        "suriname",
        "uruguay",
        "venezuela",
    ]

    americas_locations = americas_countries + [
        "new york",
        "los angeles",
        "chicago",
        "toronto",
        "vancouver",
        "montreal",
        "mexico city",
        "buenos aires",
        "são paulo",
        "sao paulo",
        "rio de janeiro",
        "lima",
        "bogotá",
        "bogota",
        "santiago",
        "caracas",
        "quito",
        "havana",
        "amazonas",
        "amazon",
        "pará",
        "para",
        "acre",
        "rondônia",
        "rondonia",
        "bolivar",
        "bolívar",
    ]

    # Europe
    europe_countries = [
        "albania",
        "andorra",
        "austria",
        "belarus",
        "belgium",
        "bosnia",
        "herzegovina",
        "bulgaria",
        "croatia",
        "czech republic",
        "czechia",
        "denmark",
        "estonia",
        "finland",
        "france",
        "germany",
        "greece",
        "hungary",
        "iceland",
        "ireland",
        "italy",
        "kosovo",
        "latvia",
        "liechtenstein",
        "lithuania",
        "luxembourg",
        "malta",
        "moldova",
        "monaco",
        "montenegro",
        "netherlands",
        "north macedonia",
        "norway",
        "poland",
        "portugal",
        "romania",
        "russia",
        "san marino",
        "serbia",
        "slovakia",
        "slovenia",
        "spain",
        "sweden",
        "switzerland",
        "ukraine",
        "united kingdom",
        "uk",
        "great britain",
        "england",
        "scotland",
        "wales",
        "northern ireland",
        "vatican",
    ]

    europe_locations = europe_countries + [
        "london",
        "paris",
        "berlin",
        "madrid",
        "rome",
        "amsterdam",
        "vienna",
        "brussels",
        "athens",
        "lisbon",
        "warsaw",
        "budapest",
        "prague",
        "stockholm",
        "moscow",
        "st petersburg",
        "kyiv",
        "kiev",
    ]

    # Pacific
    pacific_countries = [
        "australia",
        "fiji",
        "kiribati",
        "marshall islands",
        "micronesia",
        "nauru",
        "new zealand",
        "palau",
        "papua new guinea",
        "samoa",
        "solomon islands",
        "tonga",
        "tuvalu",
        "vanuatu",
    ]

    pacific_locations = pacific_countries + [
        "sydney",
        "melbourne",
        "brisbane",
        "perth",
        "auckland",
        "wellington",
        "port moresby",
        "suva",
    ]

    # Check each region with geo_loc_name
    for location in africa_locations:
        if location in location_lower:
            return "Africa", "geo"

    for location in west_asia_locations:
        if location in location_lower:
            return "West Asia", "geo"

    for location in south_asia_locations:
        if location in location_lower:
            return "South Asia", "geo"

    for location in southeast_asia_locations:
        if location in location_lower:
            return "Southeast Asia", "geo"

    for location in east_asia_locations:
        if location in location_lower:
            return "East Asia", "geo"

    for location in americas_locations:
        if location in location_lower:
            return "Americas", "geo"

    for location in europe_locations:
        if location in location_lower:
            return "Europe", "geo"

    for location in pacific_locations:
        if location in location_lower:
            return "Pacific", "geo"

    # If geo_loc_name didn't match, try isolate_name
    if isolate_name:
        isolate_lower = isolate_name.lower()

        for location in africa_locations:
            if location in isolate_lower:
                return "Africa", "isolate"

        for location in west_asia_locations:
            if location in isolate_lower:
                return "West Asia", "isolate"

        for location in south_asia_locations:
            if location in isolate_lower:
                return "South Asia", "isolate"

        for location in southeast_asia_locations:
            if location in isolate_lower:
                return "Southeast Asia", "isolate"

        for location in east_asia_locations:
            if location in isolate_lower:
                return "East Asia", "isolate"

        for location in americas_locations:
            if location in isolate_lower:
                return "Americas", "isolate"

        for location in europe_locations:
            if location in isolate_lower:
                return "Europe", "isolate"

        for location in pacific_locations:
            if location in isolate_lower:
                return "Pacific", "isolate"

    # If isolate_name didn't match, try reference_title
    if reference_title:
        reference_lower = reference_title.lower()

        for location in africa_locations:
            if location in reference_lower:
                return "Africa", "reference"

        for location in west_asia_locations:
            if location in reference_lower:
                return "West Asia", "reference"

        for location in south_asia_locations:
            if location in reference_lower:
                return "South Asia", "reference"

        for location in southeast_asia_locations:
            if location in reference_lower:
                return "Southeast Asia", "reference"

        for location in east_asia_locations:
            if location in reference_lower:
                return "East Asia", "reference"

        for location in americas_locations:
            if location in reference_lower:
                return "Americas", "reference"

        for location in europe_locations:
            if location in reference_lower:
                return "Europe", "reference"

        for location in pacific_locations:
            if location in reference_lower:
                return "Pacific", "reference"

    # If reference_title didn't match, try journal
    if journal:
        journal_lower = journal.lower()

        for location in africa_locations:
            if location in journal_lower:
                return "Africa", "journal"

        for location in west_asia_locations:
            if location in journal_lower:
                return "West Asia", "journal"

        for location in south_asia_locations:
            if location in journal_lower:
                return "South Asia", "journal"

        for location in southeast_asia_locations:
            if location in journal_lower:
                return "Southeast Asia", "journal"

        for location in east_asia_locations:
            if location in journal_lower:
                return "East Asia", "journal"

        for location in americas_locations:
            if location in journal_lower:
                return "Americas", "journal"

        for location in europe_locations:
            if location in journal_lower:
                return "Europe", "journal"

        for location in pacific_locations:
            if location in journal_lower:
                return "Pacific", "journal"

    return "Other", ""


def parse_genbank_to_csv(genbank_file, output_csv):
    """
    Parse GenBank file and write data to CSV.

    Args:
        genbank_file: Path to input GenBank file
        output_csv: Path to output CSV file
    """
    records_data = []

    # Parse GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        accession = record.id
        try:
            sequence = str(record.seq)
        except UndefinedSequenceError as e:
            sequence = ""
            print(f"Error parsing sequence: {accession}: probably contig")

        # Extract geo_loc_name and isolate from features
        geo_loc_name = ""
        isolate = ""
        for feature in record.features:
            if "geo_loc_name" in feature.qualifiers:
                geo_loc_name = feature.qualifiers["geo_loc_name"][0]
            if "isolate" in feature.qualifiers:
                isolate = feature.qualifiers["isolate"][0]

        # Extract title and journal from first reference
        reference_title = ""
        journal = ""
        if record.annotations.get("references"):
            first_ref = record.annotations["references"][0]
            if hasattr(first_ref, "title") and first_ref.title:
                reference_title = first_ref.title
            if hasattr(first_ref, "journal") and first_ref.journal:
                journal = first_ref.journal

        # Classify region - only use reference_title and journal if geo_loc_name and isolate both failed
        region, source = classify_region(
            geo_loc_name, isolate, reference_title, journal
        )

        records_data.append(
            {
                "accession": accession,
                "sequence": sequence,
                "geo_loc_name": geo_loc_name,
                "isolate": isolate,
                "reference_title": reference_title,
                "region": region,
            }
        )

    # Write to CSV
    with open(output_csv, "w", newline="") as csvfile:
        fieldnames = [
            "accession",
            "sequence",
            "geo_loc_name",
            "isolate",
            "reference_title",
            "region",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerows(records_data)

    print(f"Successfully processed {len(records_data)} records to {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_genbank_file> <output_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        parse_genbank_to_csv(input_file, output_file)
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
