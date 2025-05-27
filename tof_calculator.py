import numpy as np

IDEAL_GAS_CONSTANT = 0.0019872
KBH = 20836600000

class State():
    def __init__(self, name, value):
        self.name = name
        self.value = value
    
    def __repr__(self):
        return f'State(self.name={self.name}, self.value={self.value})'

    def __str__(self):
        return f'State {self.name}: {self.value}'
    
    def __iter__(self):
        return iter([self.name, self.value])

class CatalyticCycle():
    def __init__(self, states:list[State]):            
            self.states = states
            self.reactant = self.states[0]
            self.product = self.states[-1]
            self.ground_states = self.states[0::2]
            self.transition_states = self.states[1::2]
            self.dg = self.product.value - self.reactant.value
    
    def get_states_as_tuples(self) -> list[tuple]:
        states = []
        for state in self.states:
            states.append(tuple(state))
        
        return states
    
    def get_ground_states_as_tuples(self) -> list[tuple]:
        states = []
        for state in self.ground_states:
            states.append(tuple(state))
        
        return states
    
    def get_transition_states_as_tuples(self) -> list[tuple]:
        states = []
        for state in self.transition_states:
            states.append(tuple(state))
        
        return states
    



class ITTable():
    def __init__(self, cycle:CatalyticCycle, temperature):
        self.cycle = cycle
        self.temperature = temperature
        self.n_steps = len(cycle.states) // 2
        self.table = np.zeros((self.n_steps, self.n_steps))
        
        for row_i in range(self.n_steps):
            for col_i in range(self.n_steps):
                ts_value = self.cycle.transition_states[col_i].value
                ground_value = self.cycle.ground_states[row_i].value
                if row_i-1 < col_i:
                    dg = self.cycle.dg
                    self.table[row_i][col_i] = np.exp((ts_value - ground_value - dg) / (IDEAL_GAS_CONSTANT * self.temperature))
                else:
                    self.table[row_i][col_i] = np.exp((ts_value - ground_value) / (IDEAL_GAS_CONSTANT * self.temperature))

        self.ground_tofs = []
        self.ts_tofs = []
        for step in range(self.n_steps):
            ground_tof = self.get_row_sum(step) / self.get_total_sum()
            self.ground_tofs.append(ground_tof)

            transiton_tof = self.get_col_sum(step) / self.get_total_sum()
            self.ts_tofs.append(transiton_tof)

    def get_total_sum(self):
        return np.sum(self.table)
    
    def get_row_sum(self, index):
        return np.sum(self.table[index])
    
    def get_col_sum(self, index):
        return np.sum(self.table[:, index])

    def get_tdts(self) -> State:
        max_i = self.ts_tofs.index(max(self.ts_tofs))
        return self.cycle.transition_states[max_i]
    
    def get_tdi(self) -> State:
        max_i = self.ground_tofs.index(max(self.ground_tofs))
        return self.cycle.ground_states[max_i]

    def get_tof(self) -> float:
        return (KBH * self.temperature * 
                (np.exp(-self.cycle.dg/(IDEAL_GAS_CONSTANT * self.temperature)) - 1)) \
                / self.get_total_sum()


class TOFCalculator():
    def __init__(self, temperature:float, state_energies:list[float], state_names:list[str]=None) -> None:
        self.__temperature = temperature
        if len(state_energies) % 2 == 0:
            raise ValueError(f'state_energies must have an odd length. Current length is {len(state_energies)}')

        if state_names:
            if len(state_names) != len(state_energies):
                raise ValueError(f'state_names must have the same length as state_energies.')
            else:
                self.__state_names = state_names

        else:
            self.__state_names = []
            for i in range(len(state_energies)):
                g_or_ts = 'g' if i % 2 == 0 else 'ts'
                n = str(i // 2 + 1)
                state = g_or_ts + n
                self.__state_names.append(state)
        
        states = []
        for name, energy in zip(self.__state_names, state_energies):
            states.append(State(name, energy))
        
        self.__catalytic_cycle = CatalyticCycle(states)
        self.__it_table = ITTable(self.__catalytic_cycle, self.__temperature)

    def get_states(self) -> list[tuple]:
        return self.__catalytic_cycle.get_states_as_tuples()

    def get_tdi(self) -> tuple:
        return tuple(self.__it_table.get_tdi())

    def get_tdts(self) -> tuple:
        return tuple(self.__it_table.get_tdts())

    def get_energy_span(self) -> float:
        tdi = self.__it_table.get_tdi()
        tdts = self.__it_table.get_tdts()
        return tdts.value - tdi.value

    def get_tof(self) -> float:
        return self.__it_table.get_tof()

    def get_ground_tofs(self) -> list[tuple]:
        tofs =  self.__it_table.ground_tofs
        ground_states = [state.name for state in self.__catalytic_cycle.ground_states]
        return list(zip(ground_states, tofs))

    def get_ts_tofs(self) -> list[tuple]:
        tofs =  self.__it_table.ts_tofs
        transition_states = [state.name for state in self.__catalytic_cycle.transition_states]
        return list(zip(transition_states, tofs))


    def get_all(self) -> dict:
        return {
            'States' : self.__catalytic_cycle.get_states_as_tuples(),
            'TDI' : self.get_tdi(),
            'TDTS' : self.get_tdts(),
            'Energy Span' : self.get_energy_span(),
            'Cycle TOF' : self.get_tof(),
            'Ground State TOF Values' : self.get_ground_tofs(),
            'Transition State TOF Values' : self.get_ts_tofs()
            }


if __name__ == "__main__":
    temperature = 298.15
    state_names = ['g1', 'ts1', 'g2', 'ts2', 'g3', 'ts3', 'g4']
    state_energies = [0, 15, -7, 8, 5, 11.5, -10]
    tof_calculator = TOFCalculator(temperature=temperature,
                                   state_names=None,
                                   state_energies=state_energies)

    print(tof_calculator.get_all())


