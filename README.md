
Hidden Markov Models and Natural Language Processing


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
