# Wintermute

Analyzing and allele-calling of MPS FASTQ files. This software is part of the [DNASeqEx project](https://www.researchgate.net/project/DNASEQEX).

![European Commission Logo](https://ec.europa.eu/ec_portal/2016/images/logo/logo-splashpage.png)


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.
### Prerequisites

As this project is written Ruby, you obviously need to have a Ruby interpreter available. Please see [Install Ruby](https://www.ruby-lang.org/en/documentation/installation/) to install Ruby.

```
ruby -v
```

### Installing

First we need to install a ruby package manager called bundler. Please see [Install Bundler](https://bundler.io/) to do so.

That will give us the needed `bundle` command.

```
bundle -v
```

At the main dirctory of this project
```
bundle install
```
will install all dependencies automatically.


## Running an example

We want to analyse only the STRs of an Illumina FGx HSC FASTQ file pair (read R1 and R2 are in two different files) in the example folder of our project:

```
ls -l example
-rwx------  1 me  staff   9,8M  8 Sep  2016 D701-D501_S1_L001_R1_001.fastq.gz
-rwx------  1 me  staff   843K  8 Sep  2016 D701-D501_S1_L001_R2_001.fastq.gz
```

To call the actual program we enter

```
./wintermute.rb --only forenseq --min-reads 10 --min-call-ratio 0.1 --r2-reverse-complement --no-n-trimming example/D701-D501_S1_L001_R1_001.fastq.gz
```

The resulting `example/D501_S1_L001_CALL.xlsx` is out output file. Please note that `--r2-reverse-complement` is neccessary for Illumina FGx runs only. In reasearch mode R2 reads are already in reverse-complent order.


## Command line parameters

```
Usage: ./wintermute.rb [options] <R1 FILE>
        --[no-]calling               [Do]/[Don't do] STR/SNP calling (default: true)
        --[no-]survey                [Do]/[Don't] summarize all assigned sequences into one directory/file (default: false)
        --[no-]statistics            [Do]/[Don't] output a separate statistics file (default: false)
    -v, --[no-]verbose               Run verbosely (default: false)
    -c, --config <filename>          Configuarion file for marker and target definitions (default: "./config/generic_grch38.json")
    -k, --kits <filename>            Kit configuarion file to name markers included in kits (default: "./config/kits.json")
    -f, --[no-]force                 Overwrite result file(s) (default: false)
    -o, --output-calling FILE        STR/SNP calls output filename
    -r, --references FILE            Assign sequences to references stored in FILE
        --[no-]dynamic-q             Determine minimal Q-value dynamically based on R1/R2 (default: true)
        --[no-]n-trimming            Trim at first N (default: true)
        --append                     Append results to existing file (requires -o)
        --r1-reverse-complement      Reverse complement the R1 sequence prior use (default: false)
        --r2-reverse-complement      Reverse complement the R1 sequence prior use (default: false)
        --[no-]r2                    Do/Don't automatically determine and load R2 (default: false)
        --ignore-q                   Don't give a shit on the q-values (default: false)
        --[no-]adapter-trimming      Do/Don't trim adapter sequences (default: true)
        --[no-]primer-trimming       Do/Don't trim primer sequences (default: true)
        --require-adapter            Filter sequences that doesn't meet adapter requirements (specified in config file) (default: false)
        --only [x,y,z]               Process [x,y,z] targets/markers only. Kit names act as placeholders for the provided markers given at config/kits.json. (default: all)
        --max-distance-of-forward-sequence-from-start N
                                     Maximal distance from forward-sequence to read begin ('-1' disables this setting, default: -1)
        --[no-]reversify-targets     [Not] Reverse the configured targets to match reverse-complemnt as well (default: true)
        --no-match-forward           Sequence matching (STR/SNP) is performed only by matchng with reverse primer (default: true)
        --no-match-reverse           Sequence matching (STR/SNP) is performed only by matchng with forward primer (default: true)
        --max-n N                    Maximal number of N within matching primers and R1/R2 first 20 bases (default: 3)
        --min-q N                    Minimal Q-value for R1/R2 (overrides dynamic)
        --min-reads N                Minimal reads to consider a sequence (default: 10)
        --min-reads-ratio N          Minimal ratio of reads (relative to summarized target reads) to consider a sequence (default: 0.01)
        --min-variant-reads N        Minimal reads to consider a variant
        --min-variant-reads-ratio N  Minimal ratio of reads (relative to summarized target reads with same length) to consider a variant (default: 0.05)
        --min-call-ratio N           Minimal ratio of reads (relative to summarized target reads) to call a sequence.
```

## Configuration Files

The file `config/generic_grch38.json` contains almost al forensic relevant STR targets/markes. Please modify `config/kits.json` to fit your needs.

## Built With

* [highline](https://github.com/JEG2/highline) - A higher level command-line oriented interface.
* [rubyXL](https://github.com/weshatheleopard/rubyXL) - Ruby lib for reading/writing/modifying .xlsx and .xlsm files
* [parallel](https://github.com/grosser/parallel) - Ruby: parallel processing made simple and fast
* [pry](https://github.com/pry/pry) - An IRB alternative and runtime developer console
* [levenshtein-ffi](https://github.com/dbalatero/levenshtein-ffi) - Fast string edit distance computation, using the Damerau-Levenshtein algorithm

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Sascha Willuweit** - *Initial work*
* **Steffi Köcher** - *Testing*


## License

This project is licensed under the EUROPEAN UNION PUBLIC LICENCE v. 1.2  License - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments

The [DNASeqEx project](https://www.researchgate.net/project/DNASEQEX) has been funded with support from the European Commission (grant HOME/2014/ISFP/AG/LAWX/4000007135 under the Internal Security Funding Police programme of the European Commission-Directorate General Justice and Home Affairs).