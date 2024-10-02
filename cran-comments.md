## ReSubmission Notes 

The submission failed the pre-tests with the following notes:
package SiPhyNetwork_1.1.0.tar.gz does not pass the incoming checks automatically, please see the following pre-tests:

Windows: <https://win-builder.r-project.org/incoming_pretest/SiPhyNetwork_1.1.0_20230414_172842/Windows/00check.log>
Status: OK

Looking at the log here everything seemed fine


Debian: <https://win-builder.r-project.org/incoming_pretest/SiPhyNetwork_1.1.0_20230414_172842/Debian/00check.log>
Status: 1 NOTE

The one note had to do with a unexported function 'vcv.net'. I changed the name to vcv_net to avoid any potential method issues.


## Submission notes

The package was checked on both release and development versions of R on the latest versions of Windows, Mac, and Linux
0 errors | 0 warnings | 1 notes

There was one note from devtools:check():
checking for future file timestamps ... NOTE
  unable to verify current time

This seemed like a note that is dependent on the clock webpage: https://stackoverflow.com/questions/63613301/r-cmd-check-note-unable-to-verify-current-time


There are no known dependencies on this package

## Old Submission Notes

There was one warning from the web checks:

Version: 1.0.0
Check: whether package can be installed
Result: WARN
    Found the following significant warnings:
     code.cpp:37:10: warning: use of bitwise '|' with boolean operands [-Wbitwise-instead-of-logical]
     code.cpp:37:50: warning: use of bitwise '&' with boolean operands [-Wbitwise-instead-of-logical]
     code.cpp:46:16: warning: use of bitwise '&' with boolean operands [-Wbitwise-instead-of-logical]
Flavors: r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-fedora-clang

I made the following operands boolean instead of bitwise


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
    phylogenies (12:805)
    reticulate (12:30, 12:147, 12:280)
```
These words are correctly spelled but just uncommon

## Downstream dependencies
This is a new package and there are no downstream dependencies that I am aware of



## Old ReResubmission Notes

Thanks Benjamin for the Feedback. My notes are in bulleted points after each comment:

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking. (If you want to add a title as well please put it in
quotes: "Title")

* There are currently no references for my package. I am currently writing a manuscript but am holding off on submission until after I can get it hosted on CRAN. I will add the refrence upon submission

You have in inst/doc/introduction.R

old_pars <- par()
...
par(old_pars)

but this will produce warnings

You need:

old_pars <- par(no.readonly = TRUE)
...
par(old_pars)

* Changed this 

Please do not modify the global environment (e.g. by using <<-) in your
functions. This is not allowed by the CRAN policies.

* I do not modify the global environment by using <<-. Since the <<- is inside a function within another function, it only ever modifies the environment of the first function (foo in this case). The places where I use <<- (write.net.R and read.net.R) is essentially copy and pasted code from the 'ape' package on CRAN. My usage of the superassignment looks something like this:
```
foo<-function(){
	x<-1
	foo2<-function(y){
		x<<-x+y
	}
	
	###Do Stuff including call foo2
	foo2(3)
	foo2(4)
	###etc.
	return(x)
}
```



## Old Resubmission Notes

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


