from astroquery.splatalogue import Splatalogue


class SplatQuery:

    def __init__(self):
        self.species = species
        self.min_frequency = min_frequency
        self.max_frequency = max_frequency
        self.transition = transition
        self.show_upper_degeneracy = show_upper_degeneracy

    def lines(species,min_frequency,max_frequency,show_upper_degeneracy=True,transition=None):
        query = Splatalogue.query_lines(
                    min_frequency = min_frequency,
                    max_frequency = max_frequency,
                    chemical_name = str(" "+species+" "),
                    transition = str(transition),
                    show_upper_degeneracy = show_upper_degeneracy
                )
        return query

class LineProperties:

    def __init__(self):
        self.freq = freq

    def emissionIntensity(freq):
        P = N_u*h*omega_w*A_w