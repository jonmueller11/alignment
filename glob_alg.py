class GlobalAlignment:
    """Class to calculate the Global Alignment for to Strings.

    The class needs to sequences (as strings) and a scoring dictionary.
    After initializing the make_matrix() method needs to be called.
    Tis will create a matrix for further calculations.
    To compute the alignment the compute() method should be called.
    After that the get_alignment() method returns the alignment as a string.
    """
    def __init__(self, seq_1, seq_2, scores):
        """Initializes most attributes of an GlobalAlignment object

        :param seq_1: first sequence
        :param seq_2: second sequence
        :param scores: scores used for the alignment
                       the following key value pairs must be present:
                       'gap' : the gap penalty
                       False : value for mismatches
                       True  : value for matches
        """
        self.sequence_1: str = seq_1
        self.sequence_2: str = seq_2
        self.scores: dict = scores
        self._traceback: dict = {(0, 0): False}  # stores the traceback path
        self._seq_pos_1: list = []
        self._seq_pos_2: list = []
        self._matrix: list = []
        self._sequence_length_1: int = len(seq_1)  # the number of rows in the matrix + 1
        self._sequence_length_2: int = len(seq_2)  # the number of rows in the matrix + 1

    def make_matrix(self):
        """Creates a matrix/grid for the alignment calculation

        The method crates a matrix as a list with self._sequence_length_2 + 1
        entries. Each entry is a list of length self._sequence_length_1 + 1
        filled with NoneType objects.
        At last the first list is filled with integers cumulative summing to
        (self._sequence_length_1 + 1) * self.scores['gap']
        Every first element of each list is filled with integers cumulative summing to
        (self._sequence_length_2 + 1) * self.scores['gap']
        """
        # matrix building
        for v in range(self._sequence_length_2 + 1):
            self._matrix.append([None for _ in range(self._sequence_length_1 + 1)])
        # filling the first row and the first column
        for i in range(self._sequence_length_1 + 1):
            if i == 0:
                self._matrix[0][i] = 0
            else:
                self._matrix[0][i] = self._matrix[0][i - 1] + self.scores['gap']
        for i in range(self._sequence_length_2 + 1):
            if i == 0:
                self._matrix[i][0] = 0
            else:
                self._matrix[i][0] = self._matrix[i - 1][0] + self.scores['gap']

    def compute(self):
        """ Fills the self._matrix.

        Uses the Needleman–Wunsch algorithm[1, 2] to fill the
        self._matrix with the values.
        The tracing arrows are stored in self._traceback as key
        value pairs. Each key represents the tip of said arrow
        while the value represents the base.

        [1] https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
        [2] Saul B. Needleman, Christian D. Wunsch,
            A general method applicable to the search for similarities in the amino acid sequence of two proteins,
            Journal of Molecular Biology,
            Volume 48, Issue 3,
            1970,
            Pages 443-453,
            ISSN 0022-2836,
            https://doi.org/10.1016/0022-2836(70)90057-4.
            (https://www.sciencedirect.com/science/article/pii/0022283670900574)
        """
        for i in range(1, self._sequence_length_2 + 1):
            for j in range(1, self._sequence_length_1 + 1):
                # the calculation of the i,j value
                self._matrix[i][j] = max((
                    self._matrix[i][j - 1] + self.scores['gap'],
                    self._matrix[i - 1][j] + self.scores['gap'],
                    self._matrix[i - 1][j - 1] + self.scores[self.sequence_2[i - 1] == self.sequence_1[j - 1]]
                ))
                # adding the tracback arrow
                if self._matrix[i][j] == self._matrix[i][j - 1] + self.scores['gap']:
                    self._traceback[i, j] = (i, j-1)
                elif self._matrix[i][j] == self._matrix[i - 1][j] + self.scores['gap']:
                    self._traceback[i, j] = (i - 1, j)
                elif self._matrix[i][j] == self._matrix[i - 1][j - 1] +\
                        self.scores[self.sequence_2[i - 1] == self.sequence_1[j - 1]]:
                    self._traceback[i, j] = (i - 1, j - 1)

    def get_alignment(self):
        """Reconstructs the tracback an formats the output

        Starts with the last entry of self._matrix
        While the tracback has not reached (0, 0) each position
        is added to a list.
        After calling the self._format_alignment() method
        the method returns a formatted string, representing the alignment.
        """
        i, j = self._sequence_length_2, self._sequence_length_1
        while self._traceback[(i, j)]:
            self._seq_pos_1.append(j-1)
            self._seq_pos_2.append(i-1)
            i, j = self._traceback[(i, j)]
        self._alignment: str = self._format_alignment()
        return self._alignment

    def _format_alignment(self):
        """Returns a formatted alignment-string

        """
        # needs to reverse, so the lists are in order
        self._seq_pos_2.reverse()
        self._seq_pos_1.reverse()
        ali_1: list = []  # The Characters for the string are stored here
        ali_2: list = []  # The Characters for the string are stored here
        char_1_1: int = -1
        char_2_1: int = -1
        for char_1 in self._seq_pos_1:
            if char_1 == char_1_1:
                ali_1.append('_')  # '_' represents a gap
            else:
                ali_1.append(self.sequence_1[char_1])
            char_1_1 = char_1
        for char_2 in self._seq_pos_2:
            if char_2 == char_2_1:
                ali_2.append('_')  # '_' represents a gap
            else:
                ali_2.append(self.sequence_2[char_2])
            char_2_1 = char_2
        return ''.join(ali_1) + '\n' + ''.join(ali_2)


if __name__ == '__main__':
    sco = {'gap': 0, True: 1, False: 0}
    ali = GlobalAlignment('ACGTCATGCCCTGATGGT', 'AGTCCGGGCGCTTTAGCTAGC', sco)
    ali.make_matrix()
    ali.compute()
    print(ali.get_alignment())
