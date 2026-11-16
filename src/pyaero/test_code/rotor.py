import numpy, json

class Rotor:
    def __init__(self, parameters):
        self.parameters = parameters

    def set_design(self, parameters):
        self.parameters = parameters


with open("rotor_parameters.json", "r") as f:
    parameters = json.load(f)

rotor = Rotor(parameters)