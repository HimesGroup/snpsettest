## CRAN package review response

Fixed the issues pointed out by Gregor Seyer
* Replaced \dontrun{} with \donttest{}.
* Added reference URL enclosed in angle brackets.

## Test environments

* Arch Linux: R 4.0.4
* GitHub Actions (Ubuntu-20.04): release, devel
* GitHub Actions (Windows-lastest): release
* GitHub Actions (macOS-lastest): release
* Win-builder: devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

from Win-builder:

Possibly mis-spelled words in DESCRIPTION:
  GWAS (2:43, 11:20)
  
- It is a correctly spelled technical term (genome-wide association study)

## Downstream dependencies

There are currently no downstream dependencies for this package.
