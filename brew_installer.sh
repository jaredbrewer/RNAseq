#!/bin/sh

# This will ONLY work on Macs - macOS High Sierra (10.13) preferred, untested on other versions.

xcode-select --install
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew tap homebrew/science
brew install salmon
brew install python3
brew install R
pip3 install BioPython
pip3 install termcolor
