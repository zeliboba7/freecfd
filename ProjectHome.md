Free CFD is an open source computational fluid dynamics (CFD) code. Some of the features are detailed below but let’s first talk about why before what.

## Goals: ##

Free CFD is developed as an open source project. The idea behind the open source is not only that the code is available to everyone without any cost, but also that the code is actually understandable and usable by others. So one of the first priorities of this project is to end up with an easy-to-use, easy-to-develop code which doesn’t have too steep a learning code. That said, the project will be a sustained effort, driven by user feedback, evolving in features and interface gradually over time.

## Features: ##

### 3D Unstructured ###

Free CFD can handle arbitrary polyhedral, mixed element type 3D unstructured grids.

### Parallel ###

ParMETIS is used for domain decomposition. Open MPI is used as the message passing interface.

### All Speed ###

OK, we know that this is too general a statement but let’s say that the code can handle a Mach number of 3 as well as a Mach number of 0.001

### Density Based ###

AUSM+-up and Roe convective flux functions are currently available.

### Implicit ###

A fully impicit framework with first order, backward Euler time integration.

### Second Order Spatial Accuracy ###

Linear MUSCL reconstruction of the cell variables provide second order accuracy.