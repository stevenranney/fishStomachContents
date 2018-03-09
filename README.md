Does the recently-ingested stomach contents of fishes affect their relative weight (*W<sub>r</sub>*)? That's the question.

*** 

# NEWS

* 2018-03-09 - Manuscript sent to coauthors for review

* 2018-02-11 - All analyses are self-contained in the `stomach.R` file. The files `smallmouth.R` and `walleye.R` are gone. `source()`ing the `stomach.R` file will run all analyses and create output in a users `output` directory, provided they've created one.

* 2018-01-31 - Analyses are being maintained in the `stomach.R` file. For now, `source()`ing the `stomach.R` file will run all analyses and create the output and figures necessary for the manuscript. 
  * `smallmouth.R` and `walleye.R` are likely to go away but will remain until I've pulled salvagable code and ideas from those files.

* 2018-01-25 - Files in the `/R` dir are undergoing significant updates. Given that this code was last reviewed in June 2014, I've been doing some refactoring and combining when possible. As a result, the likelihood of `stomach.R` disappearing and the contents of `walleye.R` and `smallmouth.R` being combined is very high. `helper_functions.R` will see additions as necessary.