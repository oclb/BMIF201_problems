import numpy as np

class WrightFisherLinked:
    mutation_rate: float
    population_size: int
    selection_coefficient: float

    # This is the function that is called when you write: model = WrightFisherLinked(...)
    def __init__(self, mutation_rate: float, population_size: int, selection_coefficient: float = 0):
        self.mutation_rate = mutation_rate
        self.population_size = population_size  # N, not 2N
        self.selection_coefficient = selection_coefficient

    def scramble(self, haplotypes: np.ndarray) -> np.ndarray:
        """
        Breaks up linkage between sites by completely scrambling the carriers of each allele
        """
        for site in range(haplotypes.shape[1]):
            haplotypes[:, site] = np.random.choice(haplotypes[:, site], size=haplotypes.shape[0], replace=False)
        return haplotypes

    def next_generation(self, haplotypes: np.ndarray) -> np.ndarray:
        """
        :param haplotypes: 0/1-valued array of haplotypes where each row is a sample and each column is a site
        """

        # Effect of mutation
        mutations = np.random.choice([0, 1], size=haplotypes.shape, p=[1 - self.mutation_rate, self.mutation_rate])
        haplotypes = np.abs(haplotypes - mutations)

        # Effect of selection
        fitness = np.exp(np.sum(haplotypes, axis=1) * self.selection_coefficient)
        fitness /= np.sum(fitness)

        # Effect of drift
        samples = np.random.choice(haplotypes.shape[0], size=haplotypes.shape[0], p=fitness)

        return haplotypes[samples, :]
