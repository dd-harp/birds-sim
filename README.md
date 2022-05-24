# birds

This is a discrete time model of birds designed for simulation of West Nile virus
in bird populations, when coupled to a mosquito model.

Birds progress through an egg, fledgling, juvenile, and adult stage. Sexual
maturity is reached at the adult stage. Egg and fledgling stages are not capable
of flight. Juveniles are capable of flight, but do not mate, and consequently disperse
to find a nest. Birds capable of flight distribute their risk according to a time spent
matrix giving their home range.

A function for mating gives the time-varying (usually a yearly pulse) at which
birds seek out a mate, find a nest, and lay eggs. The model assumes these 3
events occur simultaneously, and that nest finding uses a different dispersal
matrix which will cause actual dispersal of birds.

Nestling and adult life stages can transmit disease between each other (see [Overwintering of West Nile virus in a bird community with a communal crow roost](https://www.nature.com/articles/s41598-018-24133-4)). There is no vertical
transmission. Density-dependent mortality occurs in the fledgling and juvenile stage.