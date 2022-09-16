
Thank you Victoria for initially looking over the package. Here are the comments and my attempts to address them. My notes are bulleted points after each comment.


The Description field is intended to be a (one paragraph) description of
what the package does and why it may be useful. Please add more details
about the package functionality and implemented methods in your
Description text.

* I added more text to the description.

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

* There are no references in the description

Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      biconnectedComponents.Rd: \value
      read.net.Rd: \value
      
* Added

Please write TRUE and FALSE instead of T and F.
man/network.gsa.Rd:
     network.gsa(net, ntaxa, complete = T, frac = 1, stochsampling = F)
     
* Changed this to TRUE/FALSE where I found it

Warning: Unexecutable code in man/make.trait.model.Rd:

* That code was used in an older version. I updated the code

Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
(inst/doc/introduction.R)
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)


** Reset the graphical options back to defaults

Please do not modify the global environment (e.g. by using <<-) in your
functions. This is not allowed by the CRAN policies.

** The place where I used "<<-" is a function within a function, so it should only be modifying the first functions environment. Also, the place where I see "<<-" are in 'read.net.R' and 'write.net.R' is code copied from the "ape" package that is already on CRAN. 






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
