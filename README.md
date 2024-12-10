# Tool for boosted heavy flavour jet tagger calibration

`boohft-calib` is a tool that serves for the data-MC calibration of the "boosted heavy flavour jet tagger" in CMS, based on the sfBDT coastline method.
The tool runs under the Run 2 UL condition using NanoAODv9.
It is designed for the calibration of any Xbb/Xcc type taggers composed of the branches in NanoAODv9. 

Users should specify in a data card the tagger expression, pre-defined WPs, etc., and a signal ROOT tree for extraction of the necessary signal tagger shape.
See details in the [example YMAL card](cards/example_bb_PNetXbbVsQCD.yml) for calibrating the ParticleNet XbbVsQCD score.

The introduction of the method can be found in the [BTV slides](https://indico.cern.ch/event/1120932/#23-calibration-of-ul20172018-x).
Detailed documentation is provided in [AN-21-005](https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/005) (the sfBDT method).

The calibration results and all final & intermediate plots are showcased on the webpage, automatically generated after running a routine piloted by a YAML card.
To see the example of the generated webpage, please refer to the above BTV slides.

## Run the tool

1. Run on a local cluster

First set up the environment. We recommand to use Miniconda:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ./miniconda  # for test: put the miniconda folder here
source miniconda/bin/activate
# clone the repo
git clone https://github.com/colizz/boohft-calib.git && cd boohft-calib
# install packages
conda env create -f conda_env.yml
conda activate boohft-calib
```

Run the tool in one command, e.g.,
```bash
python launcher.py cards/example_bb_PNetXbbVsQCD.yml
```

Try `python launcher.py --help` for more information on the command arguments.

Note: the tool uses 8 concurrent workers by default. On lxplus it will run by estimation 4 hrs for an entire routine. Sepcify more workers if you have more CPU resource.

2. Run on SWAN

[![Open in SWAN](https://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png)](https://cern.ch/swanserver/cgi-bin/go?projurl=https://github.com/colizz/boohft-calib.git)

Each routine will run by estimation 8 hrs on SWAN.

To run on SWAN, click the link, start a SWAN session with LCG96 Python3 stack (4 cores, 16GB), then open the `launcher_swan.ipynb` notebook and run all the blocks. This will launch a routine configured by the example card.

## Configuration card

The configuration card (e.g., the example card `cards/example_bb_PNetXbbVsQCD.yml`) defines everything for a routine. As a brief summary, users should specify
 - the type of calibration: can be `bb`, `cc`, or `qq`;
 - the year of UL condition: can be `2016APV`, `2016`, `2017`, or `2018`;
 - jet pT ranges for deriving separate SFs;
 - the tagger information, including the tagger name/expression, the span, and the custom WPs defined in the user's analysis;
 - info of a signal ROOT tree taken from the user's analysis which the tool uses for extracting the signal tagger shape.

See detailed explanation in the example card [`cards/example_bb_PNetXbbVsQCD.yml`](cards/example_bb_PNetXbbVsQCD.yml).

--------
## Update notes

v3.1.3 November 8, 2024
 - Feature: allow setting an individual fit range for the main POI
 - Feature: allow customisation of sfBDT input variables
 - Update: adapt code compatibility to the EL9 system.
 - Update: reduce the default numbers of parallel workers from 8 to 5 to prevent warnings on lxplus.

v3.1.2 July 21, 2023
 - Update: fix lumi uncertainty
 - Update: apply no JERC correction to SV mass

v3.1.1 May 25, 2023
 - Update: change the 20% frac_b/c/light variation in an overall manner (sync with mu-tagged method)
 - Update: in case of a fit failure, enlarge the autoMCStats threshold and retry
 - Feature: more text on plots to make it readable

v3.1.0 December 2, 2022
 - Feature: add new uncertainties sources
 - Feature: allow breaking down the full uncertainty list

v3.0.5 November 25, 2022
 - Feature: allow using custom sfBDT models to replace the defalt one

v3.0.4 April 19, 2022
 - Feature improved: allow expression to parse awkward-array indexing
 - Reweight binning bug fix

v3.0.3 Mar 31, 2022
 - Implement the qq calibration type

v3.0.2 Feb 5, 2022
 - Implement the year condition for 2016APV and 2016

v3.0.1 Jan 29, 2022
 - Support more command line arguments

v3.0.0 Jan 24, 2022
 - Update the method to sfBDT coastline
 - Update the framework to coffea (supports local run at present)

Previous version (till v2.1) developed in [`ParticleNet-CCTagCalib`](https://github.com/colizz/ParticleNet-CCTagCalib/)
