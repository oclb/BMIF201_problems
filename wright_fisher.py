import numpy as np


class WrightFisher:
    mutation_rate: float
    population_size: int

    # This is the function that is called when you write: model = WrightFisher(...)
    def __init__(self, mutation_rate: float, population_size: int):
        self.mutation_rate = mutation_rate
        self.population_size = population_size  # N, not 2N

    def next_generation(self, allele_frequencies: np.ndarray) -> np.ndarray:
        """
        :param allele_frequencies: array of allele frequencies in generation n
        :return: array of allele frequencies in generation n+1
        """

        expected_frequency = allele_frequencies.copy()  # avoid modifying input array

        # Effect of mutation
        # TODO: fill in the expected allele frequency in the next generation, as a function of the mutation rate
        expected_frequency =

        # Effect of drift
        # TODO: use np.random.binomial(N: integer, p: array)
        binomial_samples =

        return binomial_samples / (2 * self.population_size)

    def until_fixation(self, initial_allele_frequencies: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        :param initial_allele_frequencies: initial state
        :return: final allele frequencies (0/1-valued), time to fixation or loss
        """

        allele_frequencies = initial_allele_frequencies.copy()  # avoids modifying input array
        variants = np.where(allele_frequencies * (1-allele_frequencies) > 0)[0]  # sites not yet fixed or lost
        time_to_fixation_or_loss = np.zeros_like(allele_frequencies, dtype=int)

        generations = 0
        while len(variants) > 0:
            generations += 1
            allele_frequencies[variants] = self.next_generation(allele_frequencies[variants])
            newly_fixed_variants = (allele_frequencies[variants] * (1-allele_frequencies[variants]) == 0)
            time_to_fixation_or_loss[variants[newly_fixed_variants]] = generations
            variants = variants[~newly_fixed_variants]

        assert np.all(allele_frequencies * (1-allele_frequencies) == 0)

        return allele_frequencies, time_to_fixation_or_loss
