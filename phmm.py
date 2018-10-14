#Biopython untuk meload sekuen yang akan dialign
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

#Biopython untuk membuat HMM
from Bio.HMM import MarkovModel
from Bio.HMM import DynamicProgramming
from Bio.HMM import Trainer
from Bio.HMM import Utilities

#Profile Untuk (Diambil dari slide bioinformatik)
    # AC-SA
    # A--ST
    # ACCST
#Set Protein DNA
class DNAAlphabet(Alphabet.Alphabet):
    letters = ['A','T','S','C','']

#Set State Main, Delete, Insert
class StateAlphabet(Alphabet.Alphabet):
    letters = ['M', 'N','O','D','E','F','I','J','K','L']



#inisialisasi hmm builder
builder = MarkovModel.MarkovModelBuilder(StateAlphabet(),DNAAlphabet())

#transisi dari match state
builder.allow_transition('M','N')
builder.allow_transition('M','E')
builder.allow_transition('M','J')

builder.allow_transition('N','O')
builder.allow_transition('N','K')

builder.allow_transition('O','L')

#transisi dari insert state
builder.allow_transition('I','M')
builder.allow_transition('I','I')
builder.allow_transition('I','D')

builder.allow_transition('J','N')
builder.allow_transition('J','J')
builder.allow_transition('J','E')

builder.allow_transition('K','O')
builder.allow_transition('K','K')
builder.allow_transition('K','F')

builder.allow_transition('L','L')


#transisi dari delete state 
builder.allow_transition('D','E')
builder.allow_transition('D','J')
builder.allow_transition('D','N')

builder.allow_transition('E','F')
builder.allow_transition('E','K')
builder.allow_transition('E','O')

builder.allow_transition('F','L')


#initial probability
builder.set_initial_probabilities({'M': 1})

#probabilitas untuk transisi dari satu state ke state lain
builder.set_transition_score('M', 'N', (1/3))
builder.set_transition_score('M', 'J', (2/3))

builder.set_transition_score('J', 'J', (1/3))
builder.set_transition_score('J', 'N', (2/3))

builder.set_transition_score('N', 'O', (1))



#emission probability dari state
builder.set_emission_score('M', 'A', 1)
builder.set_emission_score('J', 'C', 1)
builder.set_emission_score('N', 'S', 1)
builder.set_emission_score('O', 'A', (1/3))
builder.set_emission_score('O', 'T', (2/3))


#building hmm
hmm = builder.get_markov_model()

#list sequence & state utk training
tseq = [Seq('ACSA',DNAAlphabet()),Seq('AST',DNAAlphabet()),Seq('ACCST',DNAAlphabet())]
tstate = [MutableSeq('MJNO',StateAlphabet()),MutableSeq('MNO',StateAlphabet()),MutableSeq('MJJNO',StateAlphabet())]

#training dengan menggunakan sequence ACSA, AST, ACCST pada tstate dan tseq
trainer = Trainer.KnownStateTrainer(hmm)
#hmm dilatih satupersatu
for i in range(len(tseq)):
    trainseq = Trainer.TrainingSequence(tseq[i],tstate[i])
    trainhmm = trainer.train([trainseq])

#Query yang akan dicocokan beserta dengan perhitungan state awalnya
seq = Seq('ACCCSA', DNAAlphabet())
state = MutableSeq('MJJJNO', StateAlphabet())

#algoritma viterbi akan mencari state yang terbaik untuk sequence
predicted_states, prob = trainhmm.viterbi(seq, StateAlphabet())

#mengeluarkan hasil probabilitas dari sequence, emission, statenya, dan predicted statenya 
print("Prediction probability: %f" % prob)
Utilities.pretty_print_prediction(seq, state, predicted_states)