# Peer-review research artifact

This repository contains the analysis code, de-identified study data, saved model objects, and generated outputs associated with the submitted manuscript.

## Contents

- `analysis.R`: analysis and figure-generation script.
- `data/raw/`: participant-level study data used by the analysis.
- `brm_median_model/`: saved Bayesian model objects and model-derived reports.
- `outputs/`: generated summaries and figures.

## Anonymisation before upload

Before upload, the authors manually prepared and reviewed a separate anonymised copy of the research archive. The following steps were applied:

- IP addresses and Prolific identifiers were replaced with `ANONYMISED`.
- Names, email addresses, external references, MTurk codes, precise timestamps, and geolocation fields were anonymised.
- Original response identifiers were replaced with randomly generated, study-specific participant identifiers using the prefixes `PILOT`, `UK`, and `USA`.
- The same participant-identifier replacements were applied inside the saved R model objects.

## Reproducing the analyses

Run `analysis.R` from the repository root. The script uses project-relative paths for the included data, models, and output directories.

## Model files

The fitted model files are included in their expected `models` directories.