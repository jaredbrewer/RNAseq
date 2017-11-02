#!/bin/sh

# This will ONLY work on Macs - macOS High Sierra (10.13) preferred, untested on other versions.

xcode-select --install
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew install salmon
pip install BioPython
pip install termcolor
