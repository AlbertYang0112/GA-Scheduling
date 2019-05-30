# GA-Homework
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Project Description
- Problem:
  Schedule the task for UAV.
- Target:
  - Minimize the length of path;
  - Minimize the execution time.
- Presumption:
  - UAV moves in a 2D plane;
  - All UAV are equal;
  - UAV is capable for all the task;
  - No limitation for the maximum number of tasks executed by a single UAV;
  - No time is needed for task execution and target is still during the task execution;
  - No consideration for the following action after the UAV finish all its task;
  - No interval between tasks.
- Constrains:
  - UAV follows the *Dubins Path*;
  - The task in sequence must be executed in order;
  - No collision between UAV.

## Our Work
- NaiveGA
  
  A simple implementation for the general GA algorithm.
  - Purpose:
    - Test the GA algorithm;
    - Explore the resolve the bottleneck of the performance;
    - An example for the GA interface;
- SchGA

  The solution for the project.
  - Task Descriptor
  - Gene Coding
  - Fitness Function
    - Dubins Path
    - Collision Detection
    - Correction for the gene violating constrains
  - Interbred
  - Mutation
  - Selection

## To Do
- bug detecting...
