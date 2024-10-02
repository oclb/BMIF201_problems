import numpy as np

class WrightFisherMigration:
    mutation_rate: float
    population_size: np.ndarray
    migration_rate: np.ndarray # square matrix whose (i,j) entry is the rate from j to i; defaults to the identity


    # This is the function that is called when you write: model = WrightFisherMigration(...)
    def __init__(self, mutation_rate: float, population_size: np.ndarray, migration_rate: np.ndarray = None):
        if migration_rate is None:
            migration_rate = np.eye(len(population_size))
        population_size = np.array(population_size)
        migration_rate = np.array(migration_rate)
        self.mutation_rate = mutation_rate
        self.population_size = population_size
        self.migration_rate = migration_rate
        self._check_migration_rate_is_valid()

    @property
    def number_of_populations(self) -> int:
        return len(self.population_size)

    def next_generation(self, allele_frequencies: np.ndarray) -> np.ndarray:
        """
        :param allele_frequencies: previous generation allele frequencies, as a number-of-populations by number-of-sites array
        :return: array of allele frequencies in the next generation
        """
        if allele_frequencies.ndim == 1:
            allele_frequencies = allele_frequencies[np.newaxis, :]

        expected_frequency = allele_frequencies.copy()  # avoid modifying input array

        # Effect of migration
        expected_frequency = self.migration_rate @ expected_frequency  # matrix multiplication

        # Effect of mutation
        expected_frequency = (expected_frequency * (1 - self.mutation_rate) +
                              (1 - expected_frequency) * self.mutation_rate)

        # Effect of drift
        result = np.zeros_like(expected_frequency)
        for i in range(self.number_of_populations):
            binomial_samples = np.random.binomial(2 * self.population_size[i], expected_frequency[i, :])
            result[i, :] = binomial_samples / (2 * self.population_size[i])

        return result

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

    def _check_migration_rate_is_valid(self) -> None:
        """
        Checks that the migration rate is a square matrix of the right size with row sums equal to 1
        """
        try:
            assert self.migration_rate.shape == (self.number_of_populations, self.number_of_populations)
        except:
            raise ValueError("Migration rate matrix must be a square matrix of size equal to the number of populations")

        row_sums = np.sum(self.migration_rate, axis=1)
        if not np.allclose(row_sums, 1):
            raise ValueError("Migration rate matrix should have row sums equal to 1")

        if not np.all(self.migration_rate >= 0):
            raise ValueError("Migration rate matrix have non-negative entries")
