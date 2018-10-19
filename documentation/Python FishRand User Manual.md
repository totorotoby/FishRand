# Python FishRand User Manual

## Inputs (Excel)
Python FishRand accepts input from a specifically formated Excel spreed sheet. The basic formatted sheet (Which can be duplicated and passed to FishRand), can be found in the FishRand subfolder:          *sheets/input/default.xlsx*
The excel spreadsheet consists of 8 different tabs, described below.
### Sample and Time input
In the time and sample tab, statistical sampling options can be set, and the time scale is defined. 
FishRand can be run in 3 different Sampling modes:

#### Deterministic mode.
In this mode, FishRand does not do Monte Carlo Simulations.  Instead, it simply uses one point value for each input parameter.  If you plan to  run FishRand with no statistical input, set both *Total number of  Uncertainty samples* and *Total number of  inner loop samples* to 1. This insures that no repeated sampling of random variables occurs.

#### Monte-Carlo Mode, but without distinguishing between variable and uncertain parameters.
To avoid distinguishing between variable and uncertain parameters, model all parameters as variable, and set *Total number of Uncertainty samples* to 1, so that the nested monte carlo simulation  enters  the  outside  uncertainty  loop  1  time  and simulates all parameters as   "variable" parameters on the inside loop.  In this mode, all statistical parameters are labeled as variable.

#### Monte-Carlo Mode with Variable and uncertain statistical input.
Both number of samples can vary and both types of variables can be defined (see "Variable Definition"). Note that since uncertain parameters are sampled from in the outer loop, if there is a low number of uncertainty samples, and a high number of variable samples final output distributions will be quite inaccurate. For good results make sure that *Total number of  Uncertainty samples* is set to at least 500.  *Total number of  Inner loop samples* can be set to 1 or more. 

#### Latin Hypercube bins

The number of Latin hypercube bins can also be set. The default is 10. With a larger number of samples fewer bins are required to give accurate output.

#### Time Input

The beginning, and end times for the simulation are entered under the sampling inputs. Beginning and end times should be entered as "MM,DD,YYYY". The time step are defined in either "Week", "Month", "Quarter", or "Day". Make sure the first letters are captialized.

#### Steady State

Lastly if the beginning and the ending times are the same, and only a single sample site is defined (if more than one is defined it will only take the first one), one can set the Steady State option to "YES" which will give a steady state solution for the model. Otherwise set Steady State to "NO". If you would like to solve for steady state with more sample sites, you can always run the number of time steps out long enough with multiple sites. 

### Parameter input formatting
In the next four input tabs, both non-statistical and statistical parameters are accepted. 

#### Adding a non-statistical parameter

To define a non-statistical parameter locate the "Entry/Distribution Parameters" pair in which you would like to input. In the "Entry" cell enter the number you would like for that parameter. Leave the "Distribution Parameters" cell completely blank (no spaces).

#### Adding a statistical parameter

Locate the "Entry/Distribution Parameters" pair. In the "Entry" cell, first define the parameter as either variable or uncertain by adding a "V" or a "U" respectively. Then define after a comma and a space, define the distribution type you would like to use. (See distribution table for distribution types).

Entry examples: 
"V, Log-Normal" 
"U, Weibull", 
"V, Beta".

Next, in the "Distribution Parameters" cell, add the corresponding distribution parameters each separated by a comma and a space (See distribution table for parameterizations).

Distribution Parameters examples: 
(For a normal distribution): 2, .5
(For a Triangle distribution): 1, 5, 4

#### Distribution Table

| Name        | Parameterization        |
| :---------- | :---------------------- |
| Normal      | μ, σ                    |
| Uniform     | beginning, length       |
| Triangle    | beginning, ending, peak |
| Log-Normal  | μ, σ *                  |
| Log-Uniform | beginning, length       |
| Beta        | α, β **                 |
| Weibull     | λ, *k*                  |

\* The Log-Normal parameters μ, and σ, are not the corresponding normal μ, and σ, but are the μ, and σ of the actual Log-Normal distribution.

\** α, and β define the beta function with pdf:  $ f(x) = \frac{\gamma(\alpha - 1) (x-1)^{b-1}}{\gamma(a)\gamma(b)} $

#### Adding Multiple Objects

In tabs 1 through 7 multiple objects (i.e. Fish, Sample Sites, Chemicals) can be defined. A new Object can be created by copy and pasting a existing object directly below itself. For instance if we wanted to have two Fish in our model, the Fish section would look like:

![Screen Shot 2018-10-04 at 1.44.10 PM](/Users/toby/Desktop/Screen Shot 2018-10-04 at 1.44.10 PM.png)

### Organism Diets

In the Org_diet tab, each Fish and Invertebrates diet must be defined. The first row of a Diet Entry specifies the organism we are created the diet for. Leave the Fraction column blank for this row. In subsequent rows a fraction of the diet can entered as a decimal point, or should be entered as 0 if this organism does not eat the other organism. The ordering of Diet Entrys must move up the food web, otherwise FishRand will crash.



### Spatial Modeling

The Last two tabs of input are for spatial modeling. If only a single Sample site is defined, or the model is in steady state these tabs can be ignored. In the Migratory_data tab at each time step, fraction of Fish present in the site can be defined as a decimal between 0 and 1. In the Sample_site_Coordinates tab the areas of each Sample Site are defined. In the first row the boundary of the entire sampling area is defined in clock wise coordinates. Below the Boundary input, a coordinate must entered for each sample site to create a thiessen polygon map. Coordinates are entered as x,y. Make sure your coordinates are within the boundary.

Lastly, Hotspots, and number of fish populations  are definied below Sample Sites. To define attraction factors set the Defintion property to either "Polygon" or "Fraction", then define its name, associated fish, attraction factor (The Density of fish relative to the outside area. For instance if a attraction factor is 2, this is equavalent to saying the density of assosiated fish in the hotspot is twice that of the outside area) and coordinates, or weights below. More coordinates or weights can be definied by adding them to the corrisponding further right columns. Make sure your hotspots are defined within the boundary.

Example:

![Screen Shot 2018-10-05 at 12.57.41 PM](/Users/toby/Desktop/Screen Shot 2018-10-05 at 12.57.41 PM.png)

Each fish must have at least one corresponding hotspot. This is a downside, but If you wish to not define a hotspot for a fish make  'fake' cooridates or weights, and set the attraction Factor to 1.



### Important Notes on Solving

**Chemical and Regional Properties**

• There are two different ways of defining Dissolved Concentration of a Chemical in Water (g/L):

​	1. You can directly define it by entry into is cell.

​	2. You can define: 

​		• Total Concentration in Water

​		•Disequilbrium factor of Disolved Organic Carbon

​		•Disequilbrium factor of Particulate Organic Carbon

​		• POC–octanol proportionality constant

​		• DOC–octanol proportionality constant

​	If 2 is used, all variables in the list must be defined.

• If any orginisims diet is partly made up of sediment, then both the Fraction of Organic Carbon Content in Sediment, and Concentration in Sediment (ng/g) must be definied.

• Some parameters in the excel spreadsheet take are labelled with deafult values that they will take on if nothing is entered.

• The model current assumes unfiltered water.

## FishRand App

The App is where you can run the excel input file, and view information about the model.

### Opening the App





### Loading Input

In the loading "Input" section of the App, use "Choose File" to find the excel input file.

 "Timesteps to Save and Display" section can be left blank if a steady state model is being ran. Otherwise the time step definied as a integer repersenting the number of timesteps since the start, can be typed in. Multiple time stops can be typed with a comma and a space between them. A time stop is the number of steps taken from the beginning date. For instance: '4, 10, 30'. Click run to run the model.

