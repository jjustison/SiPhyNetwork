## R CMD check results
The package was checked on both release and development versions of R on the latest versions of Windows, Mac, and Linux
There were no ERRORs or WARNINGs

There were two NOTEs when using `check_rhub()` to build on Windows:

```
checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

```
Possibly misspelled words in DESCRIPTION:
  Reticulate (2:37)
  reticulate (12:37)
```
These words are correctly spelled but just uncommon

## Downstream dependencies
This is a new package and there are no downstream dependencies that I am aware of
