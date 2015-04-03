AS CS440/640 Artificial Intelligence 

Hidden Markov Models and Natural Language Processing


PROBLEM DEFINITION
The goal of this assignment is to implement a basic English sentence recognizer based on Hidden Markov Models. This version has only very limited vocabulary and does not use articles or prepositions. Only the following words are included in the vocabulary supported by the HMM:

{"do" "can" "computers" "students" "movies" "develop" "play" "games" }
The HMM recognizes and parses sentences that use the above vocabulary. Our program can perform the following:
Pattern Recognition: Forward algorithm
State-Path Determination: Vertibi algorithm
Model Optimization: Baum-Welch algorithm
IMPLEMENTATION
We used Java for this project and used Eclipse as our primary IDE. The program implements the forward-backward algorithm, Vertibi algorithm and Baum-Welch algorithm. Formulas from the article by Rabiner were used to implement the algorithms. Our program can be used with command-line arguments. Here are the instructions to compile and run the program:

1. Compile: "javac HMM_System.java"
2. Run: "java HMM_System (function name) sentence.hmm (file name)"
Function names: recognize, statepath, optimize
File names: example1.obs, example2.obs, question1.obs

METHODS

Parsing Data
We began the project by parsing the HMM file given by sentence.hmm. Java's BufferedReader is used to read the data line by line. The first line of the file stores 3 variables: N (number of states), M (number of observation symbols) and T (number of length of sequence). Next the subject auxiliary predicate object words, which refer to four basic English syntactic structures, are placed into the variable SAPO. The third line contains strings which provide the vocab to be used in the observation sentences and stored in the variable list_of_vocabs. Lastly, matrix a, matrix b and matrix π are placed into separate arrays.

Pattern Recognition
For each of the observation set, the "forward" algorithm is used. In this case, the algorithm is implemented to calculate the probability of a state at a certain time given the history of evidence. The value output of the forward algorithm uses the method getAlpha, in which consists of first intializing the forward probabilities as the joint probability of state Si and the initial observation O1. 

Next the induction step solves for t>1. Using recursion, the next alpha value is calculated based on the previous alpha value and the a and b matrix values. Finally for termination, P(O|λ), which is the probability of observing the sequence O given a certain HMM, was found. Below is the formula used for the "forward" algorithm implementation. 



State Path Determination
The Viterbi algorithm is used to find the optimal state path for each observation set and reports its probabilities. In other words, the purpose of this algorithm is to find the most likely sequence of hidden states that results in a sequence of observed events. For each dataset the program recieves, the initialization, recursion and termination phases of the algorithm are implemented in the same order using the same objects defined as arrays. 

As for initialization, δ(i) is the best score (highest probability) along a single path at time t which accounts for the first t observations and ends in state S. In order to keep track of the argument for each t and j, we use array ψt(j). 

As for implementation, the recursion step that loops includes the following: loops from 1 to N through arrays δ and a to maximize the value of δ[t-1][i] * a[i][j]. The formula more in detail is shown below. 



Model Optimization
Lastly, the model is optimized using the Baum-Welch algorithm which maximizes the probability of observing a given observation set. The purpose of this optimization technique is that it incrementally produces more suitable Hidden Markov Models for a given initial model and observation sets. We used the straight-forward technique as given in the Rabiner paper for our program. However a problem with this technique is that the complexity is high because of the vast number of iterations it takes to go over the many indices in the input and intermediate values. In the future, in order to optimize the algorithm, we can create a more advanced implementation which might choose to refactor many of the calculations and repeated computation by for example, memoizing the a and b values. 

As for our initializing, given an initial set of HMM parameters in sentence.hmm and a file containing possibly multiple observation sequences, the system as a result produces as many HMM as there are input sequences. Each of the HMMs are created as a result of running just one iteration of the Baum-Welch algorithm on the corresponding input sequence. 

As for implementation, the forward-backward algorithm is used to update the parameters of λ iteratively until convergence. And in the end of an iteration, the updated values of a, b and π performs the next iteration. The result of the implementation reports back a pair of numbers: the probability of the observation sequence given the original HMM and the re-estimated HMM. Below is the formula used to optimize the HMM. 



CODE
File Name	File Description
HMM_System.java	Forward algorithm, Viterbi algorithm and Baum-Welch algorithm
sentence.hmm	HMM description
enhancement.hmm	HMM description with "Adverb" and "Present Tense"
example1.obs	Example 1 observation sequence
example2.obs	Example 2 observation sequence
question-or-statement.obs	Testing for question or statement

RESULTS

Example 1 Dataset:

1. Recognize (Forward Algorithm)
0.027
0.0288
0.0

2. Statepath (Viterbi Algorithm)
0.027 SUBJECT PREDICATE OBJECT
0.0288 SUBJECT PREDICATE OBJECT
0 

3. Optimize(Baum-Welch Algorithm)
0.027 1.0
00288 1.0
0.0 0.0

Example 2 Dataset:

1. Recognize (Forward Algorithm)
0.00058

2. Statepath (Statepath Algorithm)
0.00058 AUXILIARY SUBJECT AUXILIARY SUBJECT AUXILIARY 

3. Optimize (Baum-Welch Algorithm)
0.00058 0.037037

DISCUSSION


1. Recognize: Pattern Recognition

Why is the probability lower than we expect? 

First of all, the fact that probability of observing the output sequence is low means that it is very unlikely to observe the given observation output sequence considering the initial state possiblities, initial state transition probabilities, and initial state-output probabilities. This means that initial probabilities related to state and outputting O1 = students, O2 = play, O3 = games is very low. It could mean that the initial parameters are able to recognize only certain output sequence related with certain state sequence and state-output probabilities. This means that one hmm system can be thought of as working only for specific case of sentences.

What does this probability tell us? 

The probability tells us how likely the given observation sequence or sentence is likely to occur in the system. Therefore, we can think of it as the aptitude of the Hidden Markov Model and how well it fits the given observation. This implies that a model is usually not a general solution for all possible observation sequences. We can therefore choose one model among many models using this probability depending on which observation sequence we would like to predict.

Does the current HMM always give a reasonable answer?

No. The first and second sentences of Example 1 ("students play games" and "computers develop movies") are actual, valid sentences. But the output probability is only 2.7% and 2.88%, respectively, which are very low probability values considering that the sentences can be used in a proper English phrase. The means that the HMM is not particularly strong. The output for the third sentence ("students develop can games") has a probability of 0 which makes sense because the sentence itself holds no meaning.

What is the probability for the sentences? 
"movies do students play games" = 0.000189 = 0.0189%
"games develop play students" = 0.0 = 0.0%
2. Statepath: Model Optimization

What can we tell from the reported optimal path for syntax analysis purpose?

Given the optimal path and syntax analysis, we can determine the structure of a sentence as well as determine the purpose of the sentence, whether it is declarative or interrogative. We find the optimal path for syntax analysis purpose via state-path determination. 

Can the HMM always correctly distinguish "statement" from "question" sentence? 

Yes the HMM can always correctly distinguish a statement from a question. The reason is because in the English language a question is started with an AUXILIARY type and all we need in order to recognize if it's a question or statement is to check if the first word of the phrase is an AUXILLARY. In order to detect if the sentence is a statement, we first check if the first word is not an AUXILIARY and if we want to double check we can also check if the statement begins with a SUBJECT. I ran the ./statepath method with example phrases as shown below:

can students develop games - AUXILIARY SUBJECT PREDICATE OBJECT
students can develop games - SUBJECT AUXILIARY PREDICATE OBJECT

As seen the first statement is a question because it begins with AUXILIARY. And the second statement is a statement because it begins wtih SUBJECT.

3. Optimize: Model Optimization

Why should you not try to optimize an HMM with zero observation probability?

The main reason is because the probability of the partial observation seqeunce from t+1 to end, given state Si at time t and the model's initial parameters all become 0. Also, zero observation probability should not be used to optimize an HMM is because the Baum-Welch algorithm uses the product of the observation probabilities in the denominator of a number of intermediate values for an HMM optimization. And based on the division rule, all values cannot be divided by zero because the result will be undefined.

4. Model Enhancement

What kinds of changes will you need to make in the above HMM? Below is an example of the modified matrices a, b and π. 

For HMM models, it is very important to initialize paremeters so that we can correctly apply optimization algorithms so that we can maximize the probability of observing a given vocab sequences. As we have seen, initializing parameters with 0 makes optimization algorithm to fail to work properly. So whenver we need to initialize some probability to 0 because it is some state that is not supposed to happen in some time series, we should make it close to 0 but not 0.

Because the we are only adding on 2 more states, we will keep the original Matrix A, B and π and add values to the matrix to accomodate for "ADVERB" and "PRESENT TENSE". First we will have to change the number of states to 8 because we added new states "PRESENT TENSE" and "ADVERB" and therefore expand Matrix A. The purpose of Matrix B is the observation of the dataset. Therefore we need to add to the matrix because we have 2 new states to examine the probability of. Lastly, we expanded the π matrix by 2. However the reason why the added values are both 0.0 is because a sentence can never start with an adverb or a present tense word. Here is the result of our new sentence.hmm for the added 2 states:

Matrix A: 
0.0 0.4 0.6 0.0 0.2 0.0
0.7 0.0 0.3 0.0 0.2 0.0
0.0 0.0 0.0 0.4 0.2 0.1
0.0 0.0 0.0 0.6 0.4 0.1
0.3 0.0 0.7 0.0 0.2 0.0
0.0 0.0 0.0 0.2 0.8 0.0

Matrix B:
0.5 0.4 0.0 0.0 0.0 0.0 0.05 0.05
0.0 0.0 0.5 0.5 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.5 0.5 0.0 0.0
0.1 0.2 0.0 0.0 0.0 0.0 0.3 0.4
0.0 0.0 0.0 0.1 0.2 0.0 0.7 0.0
0.0 0.0 0.0 0.0 0.6 0.2 0.0 0.2

Matrix π:
0.6 0.3 0.1 0.0 0.0 0.0

CONCLUSION
The overall purpose of this program is to implement a simple natural language processor that is able to determine the validity of a string of English words. Although our program is able to efficiently parse the sentences, in the case of our particular implementation, it is limited to a short collection of words and probability values are lower than expected.


SOURCES
"A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition" by Lawrence Rabiner. Feb 1989.
http://www.cs.bu.edu/fac/betke/cs440/restricted/papers/rabiner.pdf

"Wikipedia: Baum Welch Algorithm"
http://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm
