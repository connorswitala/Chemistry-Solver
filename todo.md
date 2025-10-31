# TODO Plan for reorganization

Since the structure of



## Helmholtz Energy



## Gibbs Energy

I know that the first NE (number of elements) rows are reserved for elemental constraints. Here, Uj is used. Then, I have charge constraint, and finally, I have the temperature constraint if needed. So at minimum,, I will have NE correction variables (ignoring condensed phases), and at maximum I have NE + 2 correction variables. 

Clearly, the first NE rows are reserved for elemental constraints since they are always used. I will reserve row NE + 1 for charge constraint, and NE + 2 for temperature constrains. 

I can get around this by defining an integer J_SIZE and two boolean flags. One named has_ions, and one named needs_T. Then J_SIZE  = NE + has_ions + needs_T will give correct number of rows / correction variables. Then I can always say that charge constraint goes to index NE, and temperature constraint goes to J_SIZE - 1, where J_SIZE - 1 = NE if has_ions = 0.





## CESolver Library

This library should be able to solve Gibbs and Helmholtz energy. 

## Matrix Structure. 


