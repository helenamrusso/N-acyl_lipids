import re

# Polarity indices for common solvents (P' values)
polarity_index = {
    "water": 10.2,
    "methanol": 5.1,
    "ethanol": 4.3,
    "acetonitrile": 5.8,
    "isopropanol": 3.9,
    "butanol": 3.9,
    "ethyl acetate": 4.4,
    "dichloromethane": 3.1,
    "not specified": None,  # Handle unspecified cases
}

# List of solvent mixtures
all_mixtures = [
    "ethanol-water (9:1)", "ethanol-water (1:1)",
    "methanol-water (1:1)", "methanol-water (4:1)",
    "water-acetonitrile (149:1)", "methanol (100%)",
    "methanol-water (3:2)", "methanol-water (7:3)",
    "acetonitrile-water (7:3)", "acetonitrile-isopropanol-water (3:3:2)",
    "methanol-acetonitrile (3:7)", "methanol-water (9:1)",
    "acetonitrile (100%)", "acetonitrile-methanol (1:1)",
    "not specified", "dichloromethane-methanol (3:1)",
    "dichloromethane-methanol (2:1)", "ethyl acetate (100%)",
    "butanol (100%)", "ethanol-water (19:1)",
    "water (100%)", "ethanol-water (4:1)"
]

def calculate_polarity(mixture):
    match = re.match(r"([\w-]+(?:-[\w]+)*)\s*\(([\d:]+)\)", mixture)
    if match:
        solvents = match.group(1).split("-")
        ratios = list(map(int, match.group(2).split(":")))
    else:
        solvents = [mixture.replace(" (100%)", "")]
        ratios = [1]

    # Get polarity values and normalize ratios
    total = sum(ratios)
    polarity_values = [polarity_index.get(solvent, None) for solvent in solvents]

    if None in polarity_values:  # Handle unknown solvents
        return None

    # Compute weighted average polarity
    polarity = sum(p * (r / total) for p, r in zip(polarity_values, ratios))
    return polarity

# Compute and display the polarity for each mixture
polarity_results = {mix: calculate_polarity(mix) for mix in all_mixtures}

import pandas as pd
df = pd.DataFrame(list(polarity_results.items()), columns=["Mixture", "Polarity Index"])
import ace_tools as tools
tools.display_dataframe_to_user(name="Polarity Index of Mixtures", dataframe=df)
