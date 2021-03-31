import math
import random
import numpy
from numpy.distutils.fcompiler import none


def binary_search(l, x):
    left = 0
    right = len(l) - 1
    mid = right
    while left <= right:

        mid = (left + right) // 2

        if l[mid] <= x:
            left = mid + 1

        else:
            right = mid - 1

    return mid


class GeneticAlgorithm:

    def __init__(self, input_file, output_file):
        f = open(input_file, 'r')

        self.popDimension = int(f.readline())
        self.domain = tuple([int(x) for x in f.readline().split()])
        self.polynom = tuple([int(x) for x in f.readline().split()])
        self.precision = int(f.readline())
        self.setPrecision = '%.' + str(self.precision) + 'f'
        self.maxPrecision = '%.20f'
        self.crossoverProbability = float(f.readline())
        self.mutationProbability = float(f.readline())
        self.stepsNumber = int(f.readline())
        f.close()

        self.chromLength = round(math.log((self.domain[1] - self.domain[0]) * pow(10, self.precision), 2))
        self.g = open(output_file, 'w')
        self.xv = []
        self.fv = []
        self.population = []
        self.participants = []
        self.nextPopulation = []
        self.elitist_index = 0
        self.probabilities = []
        self.ruleta = []
        self.uv = []
        self.pLength = 0
        self.points = []
        self.new_xv = []
        self.new_fv = []
        self.mutants = []
        self.firstSelected = []
        self.newPopulation = []
        self.mutant_chrom = ''
        self.intervals = []

    def fitness_func(self, x):
        return self.polynom[0] * x * x + self.polynom[1] * x + self.polynom[2]

    def base2_to_10(self, chrom):
        # il aduc din baza 2 in baza 10

        chrom = int(chrom, 2)

        # si il adaptez la domeniul meu
        chrom = chrom * (self.domain[1] - self.domain[0]) / (2 ** self.chromLength - self.chromLength) \
            + self.domain[0]

        return chrom

    def elitist_element(self):
        i = (numpy.where(self.fv == numpy.max(self.fv)))[0][0]  # in liste este i,  in afisare este i + 1
        return i

    def probability(self, xi):
        s = sum(self.fv)
        return self.fv[xi] / s

    def print_population(self, chrom_list, x_list, f_list):
        for i in range(0, self.popDimension):
            # best solution?
            # discutabil: https://stackoverflow.com/questions/2189800/how-to-find-length-of-digits-in-an-integer
            sp = ' ' * (len(str(self.popDimension)) - len(str(i + 1)))
            self.g.write(str(i + 1) + sp + ': ' + str(chrom_list[i]) + '   ')

            sp = ' '
            if x_list[i] < 0:
                sp = ''

            self.g.write('x = ' + sp + self.setPrecision % (x_list[i]) + '   ')
            self.g.write('f = ' + str(f_list[i]) + '\n')

    def print_initial_population(self):
        self.g.write('Populatia initiala:' + '\n')
        self.print_population(self.population, self.xv, self.fv)

    def initial_population(self):

        for p in range(self.popDimension):
            chromosome = ''
            for i in range(self.chromLength):
                chromosome += str(random.choice([1, 0]))
            self.population.append(chromosome)

        for i in range(0, self.popDimension):
            chrom_x = self.base2_to_10(self.population[i])
            chrom_f = self.fitness_func(chrom_x)

            self.xv.append(chrom_x)
            self.fv.append(chrom_f)

    def print_selection_by_probabilities(self):
        # 1

        self.g.write('\n' + 'Probabilitatile de selectie sunt:' + '\n')

        for i in range(self.popDimension):
            sp = ' ' * (len(str(self.popDimension)) - len(str(i + 1)))
            self.g.write('cromozom ' + str(i + 1) + sp + ': ' + str(self.probabilities[i]) + '\n')

        # 2

        self.g.write('\n' + 'Intervale probabilitati de selectie:' + '\n' + '0' + ' ')

        for i in range(self.popDimension):
            self.g.write(str(self.intervals[i]) + ' ')
            if (i + 1) % 4 == 0:  # asta doar asa, ca sa fie distribuite mai uniform
                self.g.write('\n')

        # 3
        sp = ' ' * (len(str(self.popDimension)) - len(str(self.elitist_index + 1)))
        self.g.write(
            '\n' + 'selectam cromozomul ' + str(
                self.elitist_index + 1) + sp + '  prin Criteriul Elitist ' + '\n')

        # 4
        for i in range(0, self.popDimension - 1):
            sp = ' ' * (len(str(self.popDimension)) - len(str(self.ruleta[i][0] + 1)))
            self.g.write('selectam cromozomul ' + str(self.ruleta[i][0] + 1) + sp +
                         '  prin Crtietiul Ruletei u = ' + str(self.ruleta[i][1]) + '\n')

        # 5
        self.g.write('\n' + 'Dupa selectie:' + '\n')
        self.print_population(self.newPopulation, self.new_xv, self.new_fv)

    def selection_by_probabilities(self):

        # 1

        self.probabilities = []

        for i in range(self.popDimension):
            self.probabilities.append(self.probability(i))

        # 2

        self.intervals = [self.probabilities[0]]

        for i in range(1, self.popDimension):
            self.intervals.append(self.intervals[i - 1] + self.probabilities[i])

        # 3

        self.elitist_index = self.elitist_element()

        # 4

        # newPopulation incepe aici
        # n are sens sa retin doar indicii ca oricum va trebui sa fac o noua lista de cromozomi la un moment dat.

        self.newPopulation = [self.population[self.elitist_index]]
        self.new_xv = [self.xv[self.elitist_index]]
        self.new_fv = [self.fv[self.elitist_index]]
        self.ruleta = []

        while len(self.newPopulation) < self.popDimension:
            u = random.random()
            u_index = binary_search(self.intervals, u)
            self.newPopulation.append(self.population[u_index])
            self.new_xv.append(self.xv[u_index])
            self.new_fv.append(self.fv[u_index])
            self.ruleta.append((u_index, u))

    def crossover2(self, i):
        point = random.randint(0, self.chromLength - 1)

        parent1 = self.newPopulation[self.participants[i]]
        parent2 = self.newPopulation[self.participants[i + 1]]

        child1 = (parent1[:point] + parent2[point:])
        child2 = (parent2[:point] + parent1[point:])

        # replace parents with their children in new population
        self.newPopulation[self.participants[i]] = child1
        self.newPopulation[self.participants[i + 1]] = child2

        # recalculate x and f
        self.new_xv[self.participants[i]] = self.base2_to_10(child1)
        self.new_xv[self.participants[i + 1]] = self.base2_to_10(child2)
        self.new_fv[self.participants[i]] = self.fitness_func(self.xv[self.participants[i]])
        self.new_fv[self.participants[i + 1]] = self.fitness_func(self.xv[self.participants[i + 1]])

        return point

    def crossover3(self, i):
        point = random.randint(0, self.chromLength - 1)

        parent1 = self.newPopulation[self.participants[i]]
        parent2 = self.newPopulation[self.participants[i + 1]]
        parent3 = self.newPopulation[self.participants[i + 2]]

        child1 = (parent1[:point] + parent2[point:])
        child2 = (parent2[:point] + parent3[point:])
        child3 = (parent3[:point] + parent1[point:])

        # replace parents with their children in new population
        self.newPopulation[self.participants[i]] = child1
        self.newPopulation[self.participants[i + 1]] = child2
        self.newPopulation[self.participants[i + 2]] = child3

        # recalculate x and f
        self.new_xv[self.participants[i]] = self.base2_to_10(child1)
        self.new_xv[self.participants[i + 1]] = self.base2_to_10(child2)
        self.new_xv[self.participants[i + 2]] = self.base2_to_10(child3)

        self.new_fv[self.participants[i]] = self.fitness_func(self.xv[self.participants[i]])
        self.new_fv[self.participants[i + 1]] = self.fitness_func(self.xv[self.participants[i + 1]])
        self.new_fv[self.participants[i + 2]] = self.fitness_func(self.xv[self.participants[i + 2]])

        return point

    def print_crossover_stage(self):

        # 1

        self.g.write('\n' + 'Crossover probability = ' + str(self.crossoverProbability) + '\n')

        sp = ' ' * (len(str(self.popDimension + 1)) - 1)
        self.g.write(
            str(1) + sp + ': ' + self.firstSelected[0] + '  u = 0  elementul elitist nu participa la crossing over\n')

        # 2
        for i in range(1, self.popDimension):
            sp = ' ' * (len(str(self.popDimension + 1)) - len(str(i + 1)))
            self.g.write(str(i + 1) + sp + ': ' + self.firstSelected[i] + '  u = ' + str(self.uv[i]))

            if self.uv[i] <= self.crossoverProbability:
                self.g.write(' < ' + str(self.crossoverProbability) + ' Participa')

            self.g.write('\n')

        # 3

        if self.pLength <= 1:
            self.g.write('Prea putini participanti la crossing over. Nu au loc incrucisari.\n')

        if self.pLength % 2 and self.pLength > 1:
            y = 3
        else:
            y = 1

        k = 0

        for i in range(0, self.pLength - y, 2):
            self.g.write('\nIncrucisare intre cromozomul ' + str(self.participants[i] + 1) + ' si cromozomul ' +
                         str(self.participants[i + 1] + 1) + ':' + '\n')

            self.g.write(self.newPopulation[self.participants[i]] + '\n' +
                         self.newPopulation[self.participants[i + 1]] + '\npunct = ' + str(self.points[k] + 1) + '\n')
            k += 1

        if y == 3:  # mi au ramas 3 cromozomi de combinat
            self.g.write('\nIncrucisare intre cromozomii ' + str(self.participants[self.pLength - 3] + 1) + ', ' + str(
                self.participants[self.pLength - 2] + 1) + ' si ' + str(
                self.participants[self.pLength - 1] + 1) + ':' + '\n')

            self.g.write(self.newPopulation[self.participants[self.pLength - 3]] + '\n' +
                         self.newPopulation[self.participants[self.pLength - 2]] +
                         '\n' + self.newPopulation[self.participants[self.pLength - 1]] +
                         '\n' + 'punct = ' + str(self.points[k] + 1) + '\n')
            k += 1

        self.g.write('\n')

        # 4

        self.g.write('Dupa incrucisare:\n')
        self.print_population(self.newPopulation, self.new_xv, self.new_fv)

    def crossover_stage(self):

        self.firstSelected = self.newPopulation.copy()  # am nevoie de copia asta pentru afisare, atat

        # 2 probabilitatile si participantii

        self.uv = [0]

        for i in range(1, self.popDimension):
            u = random.random()
            self.uv.append(u)

            if u <= self.crossoverProbability:
                self.participants.append(i)

        # 3 shuffle + crossover

        random.shuffle(self.participants)

        self.pLength = int(len(self.participants))

        if self.pLength % 2 and self.pLength > 1:
            y = 3
        else:
            y = 1

        self.points = []

        for i in range(0, self.pLength - y, 2):
            self.points.append(self.crossover2(i))

        if y == 3:
            # mi au ramas 3 cromozomi de combinat
            self.points.append(self.crossover3(self.pLength - 3))

    def print_mutation_stage(self):

        # 1

        self.g.write('\nMutation probability ' + str(self.mutationProbability) + '\nCromozomi modificati:\n')

        for mutant in self.mutants:
            self.g.write('cromozomul ' + str(mutant[0] + 1) + ' gena ' + str(mutant[1] + 1) + '\n')

        if self.mutant_chrom == none:
            self.g.write('Niciun cromozom modificat.\n')

        # 2

        self.g.write('\nDupa mutatie:\n')
        self.print_population(self.newPopulation, self.new_xv, self.new_fv)

    def mutation_stage(self):

        # 1

        self.mutant_chrom = none
        self.mutants = []

        for i in range(1, self.popDimension):  # de la pozitia 1 pt ca nu il iau si pe elitist

            u = random.random()

            if u <= self.mutationProbability:

                p = random.randint(0, self.chromLength - 1)
                complement = str(int(not int(self.newPopulation[i][p])))

                if p == 0:
                    self.mutant_chrom = complement + \
                                        str(self.newPopulation[i][(p + 1):])

                elif p == self.chromLength - 1:
                    self.mutant_chrom = str(self.newPopulation[i][:p]) + \
                                        complement

                else:
                    self.mutant_chrom = str(self.newPopulation[i][:p]) + \
                                        complement + \
                                        str(self.newPopulation[i][(p + 1):])

                self.newPopulation[i] = self.mutant_chrom
                self.new_xv[i] = self.base2_to_10(self.mutant_chrom)
                self.new_fv[i] = self.fitness_func(self.new_xv[i])
                self.mutants.append((i, p))  # (indicele cromozomului, indicele genei)

    def evolve(self):

        self.initial_population()
        self.print_initial_population()

        self.selection_by_probabilities()
        self.print_selection_by_probabilities()

        self.crossover_stage()
        self.print_crossover_stage()

        self.mutation_stage()
        self.print_mutation_stage()

        self.g.write('\nEvolutia maximului:\n')

        for i in range(self.stepsNumber):
            self.population = self.newPopulation
            self.xv = self.new_xv
            self.fv = self.new_fv

            self.g.write(str(self.fv[self.elitist_element()]) + '\n')

            self.selection_by_probabilities()
            self.crossover_stage()
            self.mutation_stage()


def main():
    ob = GeneticAlgorithm('GeneticAlgorithm.in', 'Evolution.txt')
    ob.evolve()


main()
