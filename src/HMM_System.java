
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.OptionalDouble;
import java.text.NumberFormat;
import java.text.DecimalFormat;

public class HMM_System {
	public static void main(String[] args) {
		args = new String[4];
		args[0] = "./optimize";
		args[1] ="/Users/jungwoonshin/javaide/workspace/HMM/src/sentence.hmm";
		args[2] = "/Users/jungwoonshin/javaide/workspace/HMM/src/example1.obs";
		args[3] = "/Users/jungwoonshin/git/cs440/src/p03/output.txt";

		int i = 0;
		int N = 0;
		int M = 0;
		int T = 0;

		double[][] matrixA = null;
		double[][] matrixB = null;
		double[] matrixPI = null;

		String currentLine = null;
		String[] parsedLine = null;;
		String[] vocabList = null;
		String[] stateList = null;

		BufferedReader buffer = null;

		if (args.length > 0) {

			if ((args[0].equals("./recognize") && args.length == 3) ||
					(args[0].equals("./statepath") && args.length == 3) ||
					(args[0].equals("./optimize") && args.length == 4)) {
				try {
					File input = new File(args[1]);
					File obs = new File(args[2]);

					buffer = new BufferedReader(new FileReader(input));

					while ((currentLine = buffer.readLine()) != null) {
						parsedLine = currentLine.split(" ");

						switch(i) {
						case 0:
							N = Integer.parseInt(parsedLine[0]);
							M = Integer.parseInt(parsedLine[1]);
							T = Integer.parseInt(parsedLine[2]);

							matrixA = new double[N][N];
							matrixB = new double[N][M];
							matrixPI = new double[N];
							stateList = new String[N];
							break;
						case 1:
							stateList = currentLine.split(" ");
							break;
						case 2:
							vocabList = currentLine.split(" ");
							break;
						case 3:
							break;
						case 4:
							matrixA = extract_amatrix_Values(0, matrixA, currentLine.split(" "));
							break;
						case 5:
							matrixA = extract_amatrix_Values(1, matrixA, currentLine.split(" "));
							break;
						case 6:
							matrixA = extract_amatrix_Values(2, matrixA, currentLine.split(" "));
							break;
						case 7:
							matrixA = extract_amatrix_Values(3, matrixA, currentLine.split(" "));
							break;
						case 9:
							matrixB = extract_bmatrix_Values(0, matrixB, parsedLine);
							break;
						case 10:
							matrixB = extract_bmatrix_Values(1, matrixB, parsedLine);
							break;
						case 11:
							matrixB = extract_bmatrix_Values(2, matrixB, parsedLine);
							break;
						case 12:
							matrixB = extract_bmatrix_Values(3, matrixB, parsedLine);
							break;
						case 13:
							break;
						case 14:
							matrixPI[0] = Double.valueOf(parsedLine[0]);
							matrixPI[1] = Double.valueOf(parsedLine[1]);
							matrixPI[2] = Double.valueOf(parsedLine[2]);
							matrixPI[3] = Double.valueOf(parsedLine[3]);
							break;
						}
						i++;
					}

					// DEBUG
					//                    printMatrix(matrixA, "Matrix A Before");
					//                    printMatrix(matrixA, "Matrix B Before");
					//                    System.out.println("matrixPI before:");
					//                    System.out.println(Arrays.toString(matrixPI) + "\n");


					// Recognize Case
					if (args[0].equals("./recognize")) {
						recognize(obs, N, vocabList, matrixA, matrixB, matrixPI);

						// Statepath Case
					} else if (args[0].equals("./statepath")) {
						double[] probability = recognizeAux(obs, N, vocabList, matrixA, matrixB, matrixPI);
						statepath(obs, N, vocabList, matrixA, matrixB, matrixPI, stateList, probability);

						// Optimize Case
					} else {
						File output = new File(args[3]);


						// Optimize and extract updated matrices
						optimize(obs, vocabList, stateList, N, M, T, matrixA, matrixB, matrixPI,  output);

					}
				} catch (IOException error) {
					error.printStackTrace();
				}
			} else {
				// Error Handling
				if (args.length != 3 && (args[0].equals("./recognize"))) {
					if (args.length < 3) {
						System.out.println("Too few argument. recognize() or statepath() requires 3 args");
					} else {
						System.out.println("Too many argument. recognize() or statepath() requires 3 args");
					}
				} else if (args.length != 3 && (args[0].equals("./statepath"))) {
					if (args.length < 3) {
						System.out.println("Too few argument. recognize() or statepath() requires 3 args");
					} else {
						System.out.println("Too many argument. recognize() or statepath() requires 3 args");
					}
				} else if (args.length != 4 && args[0].equals("./optimize")) {
					if (args.length < 4) {
						System.out.println("Too few arguments. optimize() requires 4 args");
					} else {
						System.out.println("Too many arguemnts. optimize() requires 4 args");
					}
				} else {
					System.out.println("Unrecognized command: " + args[0]);
				}
			}
		}
	}


	/**************************************************** RECOGNIZE ****************************************************/


	/* 
	 * recognize() - runs the forward algorithm given inputs *.hmm and *obs and prints the probability for each dataset
	 * @param   File, int, String[], double[][], double[][] double[]
	 * @return  void
	 */
	private static void recognize(File obsFilename, int N, String[] list_of_vocabs, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix)
			throws NumberFormatException, IOException {
		BufferedReader br = new BufferedReader(new FileReader(obsFilename));
		String sCurrentLine = br.readLine();
		int dataset_length = Integer.parseInt(sCurrentLine);

		double[][] alpha;
		@SuppressWarnings("unused")
		double[][] beta;

		int numWords = 0;

		for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
			if ((sCurrentLine = br.readLine()) != null) {
				numWords = Integer.parseInt(sCurrentLine);
			}
			String words[] = br.readLine().split(" ");
			int obsIndex[] = new int[numWords];
			for (int j = 0; j < numWords; j++) {
				obsIndex[j] = getIndex(words[j], list_of_vocabs);
			}

			alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords, obsIndex);
			beta = getBeta(N, a_matrix, b_matrix, numWords, obsIndex);

			double Answer = 0.0;
			for (int state = 0; state < N; state++) {
				Answer += alpha[numWords - 1][state];
			}

			NumberFormat formatter = new DecimalFormat("#0.00000");

			System.out.println(formatter.format(Answer));
		}

		br.close();
	}


	/*
	 * recognizeAux() - A version of recognize() that returns an array of probabilities
	 * @param   File, int, String[], double[][], double[][]. double[]
	 * @return  double[]
	 */
	private static double[] recognizeAux(File obsFilename, int N, String[] list_of_vocabs, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix)
			throws NumberFormatException, IOException {
		BufferedReader br = new BufferedReader(new FileReader(obsFilename));
		String sCurrentLine = br.readLine();
		int dataset_length = Integer.parseInt(sCurrentLine);

		double[][] alpha;
		@SuppressWarnings("unused")
		double[][] beta;

		int numWords = 0;

		double[] ret = new double[dataset_length];

		for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
			if ((sCurrentLine = br.readLine()) != null) {
				numWords = Integer.parseInt(sCurrentLine);
			}
			String words[] = br.readLine().split(" ");
			int obsIndex[] = new int[numWords];
			for (int j = 0; j < numWords; j++) {
				obsIndex[j] = getIndex(words[j], list_of_vocabs);
			}
			alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords, obsIndex);
			beta = getBeta(N, a_matrix, b_matrix, numWords, obsIndex);

			double Answer = 0.0;
			for (int state = 0; state < N; state++) {
				Answer += alpha[numWords - 1][state];
			}
			ret[dataSetIndex] = Answer;
		}

		br.close();
		return ret;
	}

	/**************************************************** STATEPATH ****************************************************/


	/* 
	 * statepath() - runs the Viterbi Algorithm and prints the path with the highest probability P(O, I | lambda) and its state sequence
	 * @param   File, int, String[], double[][], double[][] double[], String[], double[]
	 * @returns void
	 */
	private static void statepath(File obsFilename, int N, String[] list_of_vocabs, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, String[] stateList, double[] probability) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(obsFilename));
		String sCurrentLine = br.readLine();
		int dataset_length = Integer.parseInt(sCurrentLine);

		@SuppressWarnings("unused")
		double[][] alpha;
		@SuppressWarnings("unused")
		double[][] beta;

		int numWords = 0;

		for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
			if ((sCurrentLine = br.readLine()) != null) {
				numWords = Integer.parseInt(sCurrentLine);
			}
			String words[] = br.readLine().split(" ");
			int obsIndex[] = new int[numWords];
			for (int j = 0; j < numWords; j++) {
				obsIndex[j] = getIndex(words[j], list_of_vocabs);
			}
			alpha = getAlpha(N, a_matrix, b_matrix, pi_matrix, numWords, obsIndex);
			beta = getBeta(N, a_matrix, b_matrix, numWords, obsIndex);
			runViterabi(N, list_of_vocabs, a_matrix, b_matrix, pi_matrix, numWords, obsIndex, stateList, dataSetIndex, probability);
		}
		br.close();
	}

	/* 
	 * runViterabi - Viterabi Algorithm
	 * @param   int, Stringp[, double[][], double[][], double[], int, int[], Stringp[, int, double[]
	 * returns  void
	 */
	public static void runViterabi(int N, String[] list_of_vocabs,
			double[][] a_matrix, double[][] b_matrix, double[] pi_matrix,
			int numWords, int[] obsIndex, String[] stateList, int index, double[] probability) {
		double[][] sigma = new double[numWords][N];
		double[][] chai = new double[numWords][N];
		double[] a_ij_times_sigma = new double[N];

		// Step 1 & 2: Recursion
		sigma = getSigma(N, a_matrix, b_matrix, pi_matrix, sigma, numWords, obsIndex, a_ij_times_sigma);
		chai = getChai(N, a_matrix, sigma, chai, numWords, a_ij_times_sigma);

		// Step 3: termination step
		int q_star[] = new int[numWords];
		q_star = getQStar(sigma, chai, numWords, q_star);

		NumberFormat formatter = new DecimalFormat("#0.00000");
		System.out.print(formatter.format(probability[index]) + " ");

		for (int i = 0; i < q_star.length; i++) {
			System.out.print(stateList[q_star[i]] + " ");
		}
		System.out.println();
	}


	/**************************************************** OPTIMIZE ****************************************************/


	/* 
	 * optimize() - runs the Baum-Welch Algorithm to optimize the input and returns a list of updated matrices
	 * @param   File, String[], int, int, double[][], double[][], double[]
	 * @returns ArrayList<double[][]
	 */
	private static void optimize(File obsFilename, String[] vocabList, String[] stateList, int N, int M, int T, double[][] matrixA, double[][] matrixB, double[] matrixPI, File output)
			throws NumberFormatException, IOException {
		BufferedReader br = new BufferedReader(new FileReader(obsFilename));
		String sCurrentLine = br.readLine();
		int dataset_length = Integer.parseInt(sCurrentLine);

		double[][] alpha;
		double[][] beta;
		int numWords = 0;
		//		System.out.println("dataset_length: " + dataset_length);
		for (int dataSetIndex = 0; dataSetIndex < dataset_length; dataSetIndex++) {
			if ((sCurrentLine = br.readLine()) != null) {
				numWords = Integer.parseInt(sCurrentLine);
			}
			String words[] = br.readLine().split(" ");
			int obsIndex[] = new int[numWords];
			for (int j = 0; j < numWords; j++) {
				obsIndex[j] = getIndex(words[j], vocabList);
			}
			alpha = getAlpha(N, matrixA, matrixB, matrixPI, numWords, obsIndex);
			beta = getBeta(N, matrixA, matrixB, numWords, obsIndex);

			runBaumWelchFix(obsFilename, vocabList, stateList, N, M, T, matrixA, matrixB, matrixPI, alpha, beta, numWords, obsIndex,  output);

		}

		br.close();
	}

	/* Macro Definition */
	private static final double DENOMINATOR_CONSTANT = 10e-70;

	public static void runBaumWelchFix(File obsFilename, String[] vocabList, String[] stateList, int N, int M, int T, double[][] a_matrix, double[][] b_matrix,
			double[] pi_matrix, double[][] alpha, double[][] beta, int numWords, int[] obsIndex, File output) throws NumberFormatException, IOException {
		double[][][] xi;
		double[][] gamma;

		xi = getXI(N, a_matrix, b_matrix, alpha, beta, numWords, obsIndex);
		gamma = getGamma(N, numWords, alpha, beta, xi);

		for (int i = 0; i < N; i++) {
			pi_matrix[i] = gamma[0][i];
		}

		double[][] trained_a_matrix = new double[N][N];

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				for (int t = 0; t < numWords - 1; t++) {
					trained_a_matrix[i][j] += xi[t][i][j];
				}
			}
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double denom = 0.0;

				for (int t = 0; t < numWords - 1; t++) {
					denom += gamma[t][i];
				}

				if (denom == 0.) {
					denom = DENOMINATOR_CONSTANT;
				}
				trained_a_matrix[i][j] /= denom;
			}
		}
		a_matrix = trained_a_matrix;

		double[][] trained_b_matrix = new double[N][M];

		for (int i = 0; i < N; i++) {
			for (int m = 0; m < M; m++) {
				for (int t = 0; t < numWords; t++) {
					if (obsIndex[t] == m) {
						trained_b_matrix[i][m] += gamma[t][i];
					}
				}
			}
		}
		for (int i = 0; i < N; i++) {
			for (int m = 0; m < M; m++) {
				double denom = 0.0;
				for (int t = 0; t < numWords; t++) {
					denom += gamma[t][i];
				}
				if (denom == 0.0) {
					denom = DENOMINATOR_CONSTANT;
				}
				trained_b_matrix[i][m] /= denom;
			}
		}
		b_matrix = trained_b_matrix;
		double unoptimized = 0.0;
		for (int state = 0; state < N; state++) {
			unoptimized += alpha[numWords - 1][state];
		}

		alpha = getAlpha(N, trained_a_matrix, trained_b_matrix, pi_matrix, numWords, obsIndex);
		beta = getBeta(N, trained_a_matrix, trained_b_matrix, numWords, obsIndex);

		double optimized = 0.0;
		for (int state = 0; state < N; state++) {
			optimized += alpha[numWords - 1][state];
		}
		System.out.println(unoptimized + " " +optimized);
		//		printMatrix(a_matrix, "Matrix A after optimization");
		//		printMatrix(b_matrix, "Matrix B after optimization");

		//		recognize(obsFilename, N,vocabList,a_matrix,b_matrix,pi_matrix);



		// Write optimized hmm to output file
		BufferedWriter fw = new BufferedWriter(new FileWriter(output));

		fw.write(N + " " + M + " " + T);
		fw.newLine();

		fw.write(Arrays.toString(stateList));
		fw.newLine();

		fw.write(Arrays.toString(vocabList));
		fw.newLine();

		// Write Updated Matrix A
		fw.write("a:");
		fw.newLine();
		writeMatrix(fw, a_matrix);

		// Write Updated Matrix B
		fw.write("b:");
		fw.newLine();
		writeMatrix(fw, b_matrix);

		// Write Updated Matrix PI
		fw.write("pi:");
		fw.newLine();
		for(int p = 0; p < pi_matrix.length; p++) {
			fw.write(pi_matrix[p] + " ");
		}
		fw.newLine();
		fw.flush();
		fw.close();


	}


	//    public static void runBaumWelch(int N, int M, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, double[][] alpha, double[][] beta,
	//                                    int numWords, int[] obsIndex) {
	//        double[][][] xi = getXI(N, a_matrix, b_matrix, alpha, beta, numWords, obsIndex);
	//        double[][] gamma = getGamma(N, numWords, alpha, beta);
	//        
	//        for (int i = 0; i < N; i++) {
	//            pi_matrix[i] = gamma[0][i];
	//        }
	//        
	//        double[][] trained_a_matrix = new double[N][N];
	//        
	//        for (int i = 0; i < N; i++) {
	//            for (int j = 0; j < N; j++) {
	//                for (int t = 0; t < numWords - 1; t++) {
	//                    trained_a_matrix[i][j] += xi[t][i][j];
	//                }
	//            }
	//        }
	//        
	//        for (int i = 0; i < N; i++) {
	//            for (int j = 0; j < N; j++) {
	//                double denom = 0.0;
	//                
	//                for (int t = 0; t < numWords - 1; t++) {
	//                    denom += gamma[t][i];
	//                }
	//                
	//                if (denom == 0.) {
	//                    denom = DENOMINATOR_CONSTANT;
	//                }
	//                trained_a_matrix[i][j] /= denom;
	//            }
	//        }
	//        
	//        double[][] trained_b_matrix = new double[N][M];
	//        
	//        for (int i = 0; i < N; i++) {
	//            for (int m = 0; m < M; m++) {
	//                for (int t = 0; t < numWords; t++) {
	//                    if (obsIndex[t] == m) {
	//                        trained_b_matrix[i][m] += gamma[t][i];
	//                    }
	//                }
	//            }
	//        }
	//        for (int i = 0; i < N; i++) {
	//            for (int m = 0; m < M; m++) {
	//                double denom = 0.0;
	//                for (int t = 0; t < numWords; t++) {
	//                    denom += gamma[t][i];
	//                }
	//                if (denom == 0.0) {
	//                    denom = DENOMINATOR_CONSTANT;
	//                }
	//                trained_b_matrix[i][m] /= denom;
	//            }
	//        }
	//
	//        alpha = getAlpha(N, trained_a_matrix, trained_b_matrix, pi_matrix, numWords, obsIndex);
	//        beta = getBeta(N, trained_a_matrix, trained_b_matrix, numWords, obsIndex);
	//    }


	/**************************************************** HELPER FUNCTIONS ****************************************************/


	/*
	 * getAlpha() -
	 * @param   int, double[][], double[][], double[], int, int[]
	 * @return  double[][]
	 */
	public static double[][] getAlpha(int N, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, int numWords, int[] obsIndex) {
		double[][] alpha = new double[numWords][N];
		// Step : t = 1
		for (int state = 0; state < N; state++) {
			alpha[0][state] = pi_matrix[state] * b_matrix[state][obsIndex[0]];
		}

		// Step : t = 2
		for (int t = 1; t < numWords; t++) {
			for (int l = 0; l < N; l++) {
				for (int k = 0; k < N; k++) {
					alpha[t][l] += alpha[t - 1][k] * a_matrix[k][l];
				}
				alpha[t][l] *= b_matrix[l][obsIndex[t]];

			}
		}
		return alpha;
	}

	/*
	 * getBeta() -
	 * @param   int, double[][], double[][], int, int[]
	 * @return  double[][]
	 */
	public static double[][] getBeta(int N, double[][] a_matrix, double[][] b_matrix, int numWords, int[] obsIndex) {
		/*
		double[][] beta = new double[numWords][N];
		for (int t = numWords; t > 0; t--) {
			for (int i = 0; i < N; i++) {
				if (t == numWords) {
					beta[numWords - 1][i] = 1;
				} else {
					for (int j = 0; j < N; j++) {
						beta[t - 1][i] += a_matrix[i][j] * b_matrix[j][obsIndex[t]] * beta[t][j];
					}
				}

			}
		}
		return beta;
		 */
		double[][] beta;
		beta = new double[numWords][N];

		// Initialization
		for (int i = 0; i < N; i++) {
			beta[numWords - 1][i] = 1.0;
		}

		// Induction
		for (int t = numWords - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {
				double sum = 0.0;
				for (int j = 0; j < N; j++) {
					sum += a_matrix[i][j] * b_matrix[j][obsIndex[t + 1]] * beta[t + 1][j];
				}
				beta[t][i] = sum;
			}
		}
		return beta;
	}

	/*
	 * getChai() -
	 * @param   int, double[][], double[][], double[][], int, double[]
	 * @return  double[][]
	 */
	public static double[][] getChai(int N, double[][] a_matrix, double[][] sigma, double[][] chai, int numWords, double[] a_ij_times_simga) {

		for (int t = 1; t < numWords; t++) {
			for (int j = 0; j < N; j++) {
				for (int i = 1; i < N; i++) {
					a_ij_times_simga[i] = sigma[t - 1][i] * a_matrix[i][j];
				}
				double highest_value = Arrays.stream(a_ij_times_simga).max().getAsDouble();

				int highest_index = 0;
				for (double s : a_ij_times_simga) {
					if (s == highest_value) {
						break;
					}
					highest_index++;
				}
				chai[t][j] = highest_index;
			}
		}
		return chai;
	}

	/*
	 * getGamma() -
	 * @param   int, int, double[][], double[][]
	 * @return  double[]
	 */
	public static double[][] getGamma(int N, int numWords, double[][] alpha, double[][] beta, double[][][] xi) {

		double[][] gamma = new double[numWords][N];

		for (int t = 0; t < numWords; t++) {
			for (int i = 0; i < N; i++) {
				gamma[t][i] = gamma(i, t, N, alpha, beta);
			}
		}
		return gamma; 
	}

	/*
	 * gamma()
	 * @param   int, int, int, double[][], double[][]
	 * @return  double
	 */
	public static double gamma(int i, int t, int N, double[][] fwd, double[][] bwd) {

		double num = fwd[t][i] * bwd[t][i];
		double denom = 0;

		for (int j = 0; j < N; j++) {
			denom += fwd[t][j] * bwd[t][j];
		}

		return divide(num, denom);
	}

	/*
	 * divide()
	 * @param   double, double
	 * @return  double
	 */
	public static double divide(double n, double d) {
		if (n == 0) {
			return 0;
		} else {
			return n / d;
		}
	}

	/*
	 * getIndex() - finds the index that matches the vocabulary list given a search string
	 * @param   String, String
	 * @return  int
	 */
	public static int getIndex(String word, String vocab[]) {
		int i = 0;
		while (i < vocab.length && !vocab[i].equalsIgnoreCase(word)) {
			i++;
		}
		return i;
	}

	/*
	 * getQStar() -
	 * @param   double[][], double[][], int, int[]
	 * @return  int[]
	 */
	public static int[] getQStar(double[][] sigma, double[][] chai, int numWords, int[] q_star) {
		double p_star = Arrays.stream(sigma[numWords - 1]).max().getAsDouble();

		for (double e : sigma[numWords - 1]) {
			if (e == p_star) {
				break;
			}
			q_star[numWords - 1]++;
		}

		// Step 4:  Path
		for (int i = numWords - 2; i >= 0; i--) {
			q_star[i] = (int) chai[i + 1][(int) q_star[i + 1]];
		}
		return q_star;
	}

	/*
	 * getSigma() -
	 * @param   int, double[][], double[][], double[], double[][], int, int[], double[]
	 * @return  double[][]
	 */
	public static double[][] getSigma(int N, double[][] a_matrix, double[][] b_matrix, double[] pi_matrix, double[][] sigma, int numWords,
			int[] obsIndex, double[] a_ij_times_simga) {
		// Initialization Step: t = 1
		for (int state = 0; state < N; state++) {
			sigma[0][state] = pi_matrix[state] * b_matrix[state][obsIndex[0]];
		}

		// Recursion Step, t = 2
		for (int t = 1; t < numWords; t++) {
			for (int j = 0; j < N; j++) {
				for (int state = 0; state < N; state++) {
					a_ij_times_simga[state] = sigma[t - 1][state] * a_matrix[state][j];
				}
				OptionalDouble highest = Arrays.stream(a_ij_times_simga).max();
				sigma[t][j] = highest.getAsDouble() * b_matrix[j][obsIndex[t]];
			}
		}
		return sigma;
	}

	/*
	 * getXI() -
	 * @param   int, double[][], double[][], double[][], double[][], int, int []
	 * @return  double[][][]
	 */
	public static double[][][] getXI(int N, double[][] a_matrix, double[][] b_matrix, double[][] alpha, double[][] beta, int numWords, int[] obsIndex) {
		double[][][] xi = new double[numWords][N][N];
		for (int t = 0; t < numWords - 1; t++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					xi[t][i][j] = getXi_IJ_value(t, i, j, N, a_matrix, b_matrix, alpha, beta, xi, numWords, obsIndex);
				}
			}
		}
		return xi;
	}

	/* 
	 * getXI_IJ_Value() -
	 * @param   int, int, int, int, double[][], double[][], double[][], double[][], double[][][], int, int[]
	 * @return  double
	 */
	public static double getXi_IJ_value(int t, int state1, int state2, int N, double[][] a_matrix, double[][] b_matrix, double[][] alpha,
			double[][] beta, double[][][] xi, int numWords, int[] obsIndex) {
		xi[t][state1][state2] = alpha[t][state1] * a_matrix[state1][state2] * b_matrix[state2][obsIndex[t + 1]] * beta[t + 1][state2];

		// Wiki version
		double denominator = 0.;
		for (int k = 0; k < N; k++) {
			denominator += alpha[t][k] * beta[t][k];
		}

		if (denominator == 0.) {
			denominator = DENOMINATOR_CONSTANT;
		}
		xi[t][state1][state2] /= denominator;

		return xi[t][state1][state2];
	}

	/*
	 * extract_amatrix_Values() - Reads matrix A and returns N*N array
	 * @param int, double[][], String[]
	 * @return double[][]
	 */
	public static double[][] extract_amatrix_Values(int i, double[][] a_matrix, String[] parseSpaceAvalues) {
		for (int k = 0; k < a_matrix[0].length; k++) {
			a_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
		}
		return a_matrix;
	}

	/*
	 * extract_bmatrix_Values() - Read matrix B and returns N*M array. N is number of states
	 * @param int, double[][] String[]
	 * @return double[][]
	 */
	public static double[][] extract_bmatrix_Values(int i, double[][] b_matrix, String[] parseSpaceAvalues) {
		for (int k = 0; k < b_matrix[0].length; k++) {
			b_matrix[i][k] = Double.valueOf(parseSpaceAvalues[k]);
		}
		return b_matrix;
	}

	/* 
	 * printMatrix() - prints the contents of a 2D matrix (Debugging Helper)
	 * @param   double[][], String
	 * @return  void
	 */
	public static void printMatrix(double[][] matrix, String matrixName) {
		System.out.println(matrixName);
		for (int row = 0; row < matrix.length; row++) {
			for (int col = 0; col < matrix[0].length; col++) {
				System.out.print(matrix[row][col] + " ");

				if (col == matrix[0].length - 1) {
					System.out.println();
				}
			}
		}
		System.out.println();
	}

	/* 
	 * writeMatrix() - writes a matrix to a specified BufferdWriter
	 * @param   BufferWriter, double[][]
	 * @return  void
	 */
	public static void writeMatrix(BufferedWriter fw, double[][] matrix) {
		for (int row = 0; row < matrix.length; row++) {
			for (int col = 0; col < matrix[0].length; col++) {
				try {
					fw.write(matrix[row][col] + " ");

					if (col == matrix[0].length - 1) {
						fw.newLine();
					}
				} catch (IOException error) {
					error.printStackTrace();
				}
			}
		}
	}
}
