# TurnOverFrequencyCalculator
A python tool for calculating turn over frequency, energy span, TOF-determining transition state and TOF-determining intermediate for homogenous catalytic cycles.  

Based on the fortran program presented in _Theoretical Analysis of the Catalytic Cycle of a Nickel Cross-Coupling Process: Application of the Energetic Span Model_ by Kozuch et al. https://doi.org/10.1021/om800772g  
## Dependencies   
- Numpy
```python
pip install numpy
```
# Documentation
## TOFCalculator  
Provides access to TOF, energy span, TDTS, TDI, and TOF values for transition states and ground states.  
### Import and Instantiation  
#### \_\_init\_\_(state_energies='', state_names=None, temperature=298)  
>PARAMETERS:  
>>- state_energies _list[float]_: list of energies of cycle states in order. Length must be odd numbered.  
>>- state_names _list[str]_: list of names of cycle states in order. Names will be generated if not provided.  
>>- temperature _float_: Temperature of reaction.  
```python  
from tof_calculator import TOFCalculator

temperature = 298
state_names = ['g1', 'ts1', 'g2', 'ts2', 'g3', 'ts3', 'g4']
state_energies = [0, 15, -7, 8, 5, 11.5, -10]

tof_calculator = TOFCalculator(temperature=temperature,
                               state_names=state_names,
                               state_energies=state_energies)
```
### Methods  
**get_states() -> _list[tuple]_**  
>returns a list of (name, energy) pairs for each state.  

**get_tdts() -> _tuple_**  
>returns name and energy of TOF determining transition state.  

**get_tdi() -> _tuple_**  
>returns name and energy of TOF determining ground state.  

**get_energy_span() -> _float_**  
>returns the difference between TDTS and TDI energies.  

**get_tof() -> _float_**  
>returns the turn over frequency of the catalytic cycle.  

**get_ground_tofs() -> _list[tuple]_**  
>returns pairs of (name, degree of TOF control) for each ground state.    

**get_ts_tofs() -> _list[tuple]_**  
>returns pairs of (name, degree of TOF control) for each transition state.  

**get_all() -> _dict_**  
>returns all previous method outputs in dictionary format.  
