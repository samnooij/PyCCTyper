#!/usr/bin/env python3

from pathlib import Path

original_profiles = list(Path("data", "Profiles").glob("*.hmm"))

output_dir = Path("pyhmmer_profiles")
output_dir.mkdir(exist_ok=True)


for original_profile in original_profiles:
    profile_name = original_profile.stem
    output_file = output_dir / original_profile.name
    with open(output_file, "w") as outfile:
        with open(original_profile, "r") as infile:
            for line in infile:
                if line.startswith("NAME"):
                    line = f"NAME  {profile_name}\n"
                outfile.write(line)
