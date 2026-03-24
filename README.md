# NESS - Name TBD

## Overview
Liquid rocket engine sizing and regenerative cooling analysis and design. It can create performance maps across throttle conditions for engine thermals and performance.

## Features
- Design an engine contour and get performance values (using RocketIsp)
- Design the thermals/fluid of regenerative cooling circuit

## Future Features
- Basic film cooling model 
- Variable area regen coolant channels
- Stress analysis on liner hot wall
- Perform parametric/sensitivity studies on channel geometry to analyze liner thermal/fluids behavior
- Two-pass cooling

RocketIsp is used for all performance calculations. See those docs here - https://rocketisp.readthedocs.io/en/latest/ 

## Structure
Explain folders briefly

## Installation on Windows 10
Many external libraries are used that must be downloaded.

## Verify Installation
To verify that all packages were installed correctly, run "" in your terminal.

DO NOT MODIFY any files yet as you might get different results. This file will size a single thrust chamber and design a regenerative cooling circuit for it. The file should take around X minutes to run. 

git clone ...
pip install -r requirements.txt

## Usage

## Documentation
The math used in this software is documented here - https://www.notion.so/Regen-Code-Overview-2a227990d78e80a1baedc85dda9b06e2?source=copy_link  

This Notion page contains detailed documentation on all heat transfer models and explains the equations in-depth. The code can be run without reading this page but please read through this to understand the functionality. Remember that this is an analysis tool and will not be useful if the analyst doesn't understand the fundamentals.

## Code Validation
Thermal/fluids/structural/performance analysis will be compared to test data. The intent of this tool is to be trusted and verified by analysts to minimize the repetitive of making new codes. Most of the time, users assume the code is a blackbox and make their own. This code will be very thoroughly documented such that users can feel confident in the results.

See code validation here: [insert]