import java.util.PriorityQueue;
import java.util.Random;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Implementation of a game board for Sliding 8-Puzzles to track states.
 * Change heuristic by changing the constant parameter field HEURISTIC.
 * 
 * @author James Ni
 * @author Arthur Chen
 */
public class SlidingGameState implements Comparable<SlidingGameState>{

	private final int DIMENSION = 3; 	//Dimension of puzzle
	private final int HEURISTIC = 3;	//Value of 1, 2, or 3 to choose heuristic. Defaults to heuristic 3
	private final boolean DEPTH = true; //Depth based equal expected cost resolution	
	private final boolean HRESOLUTION = false;	//Heuristic based resolution. If DEPTH==true; will use DEPTH
	private final static int TRIALS = 200; 	//Number of data points to take. For consistency

	private int[] state;
	private int zero;
	private int h1;
	private int h2;
	private int h3;
	private int depth;

	/**
	 * Main function to test individual puzzles
	 * @param args Input at compile time if desired
	 */
	public static void main(String[] args){

		//		int[] t0 = {0,2,1,7,4,5,6,3,8};
		//		int[] t = {2,7,0,5,4,3,8,1,6};
		//		int[] tt = {4,3,6,8,0,7,5,2,1};
		//		int[] ttt = {0,2,1,5,4,3,6,7,8};
		//		int[] tttt = {0,1,2,5,3,4,6,7,8};
		//		SlidingGameState test2 = new SlidingGameState(ttt,0,0,0);
		//		System.out.println(test2.h2);
		//		System.out.println(test2.h3);
		//		System.out.println(new SlidingGameState(24));



		double start = System.nanoTime(); //Run time testing

		//To have a consistent number of runs per path length
		ArrayList<PriorityQueue<Integer> > avg = new ArrayList<PriorityQueue<Integer> >(); 
		for(int i = 0; i<12; i++){
			avg.add(new PriorityQueue<Integer>());	//Initialie each PQ
		}

		//		int rounds = 0;
		System.out.println("Trials run per depth: " + TRIALS);
		for (int i = 24; i > 0 ; i-=2){
			PriorityQueue<Integer> pq = avg.get(i/2-1); //For simpler reference

			while(pq.size() < TRIALS){ //Run trials until the current, largest path size is not filled.
				SlidingGameState test = new SlidingGameState(i); 	//Generate a random initial state
				int[] results = test.solveBoard();
				int index = results[0]/2-1;							//Index the result to the correct PQ
				if(index >= 0 && avg.get(index).size()<TRIALS){
					avg.get(index).add(results[1]);					//Add the number of nodes expanded to be summed
				}
				//				rounds++;


			}	
			//			System.out.println(rounds);
		}

		//Sums the nodes for each path length to find an average.
		for (int i=0; i<12; i++){
			int average = 0;
			PriorityQueue<Integer> pq = avg.get(i);
			//			System.out.println(pq.size());
			while(!pq.isEmpty()){
				average+=pq.poll();
			}
			System.out.println("For solution depth equals to " + ((i+1)*2) + " there was an average of " + 
					average/TRIALS + " nodes expanded");
		}

		//Prints total execution time
		System.out.println("Execution Time (s): " + (System.nanoTime()-start)/1000000000);


	}

	/**
	 * Game State creator with inputed parameters. Used for children state generation.
	 * 
	 * @param s The board state
	 * @param depth The depth of the node
	 * @param optimal The maximum desired path to the node
	 * @param zero Current position of the empty space.
	 */
	private SlidingGameState(int [] s, int depth, int zero){
		state = s;
		h1 = getHeuristic1();
		h2 = getHeuristic2();
		h3 = getHeuristic3();
		this.depth = depth;
		this.zero = zero;
	}


	/**
	 * Initial game state creator, with ideally the inputed amount of tile shifts to reach the
	 * goal state. Can be lower as some combinations of moves will undo previous shifts.
	 * 
	 * @param moves The maximum number of moves to reach the goal state.
	 */
	public SlidingGameState(int moves){
		if(moves<0){
			throw new IllegalArgumentException("Number of moves must be non-negative.");
		}
		state = new int[DIMENSION*DIMENSION];
		Random rand = new Random(System.nanoTime());
		for(int i=0; i<DIMENSION*DIMENSION; i++){
			state[i]=i;
		}
		int shifts = 0; //Number of moves done
		int previous = 0; //Avoids undoing the last move to get more consistent optimal move results
		//Pseudo-randomly generate movement options while avoiding going out of bounds 
		// And avoiding undoing the previous movement. 
		while(shifts < moves){ 
			double direction = rand.nextDouble(); 
			//If it is not on the left edge and doesn't undo
			if(direction < .25 && zero%DIMENSION != 0 && previous!=1){
				move(-1);
				previous = -1;
				shifts++;
			}
			//If it is not on the right edge and doesn't undo
			else if(direction < .5 && zero%DIMENSION != DIMENSION-1 && previous!=-1){
				move(1);
				previous = 1;
				shifts++;
			}
			//If it is not on the top edge and doesn't undo
			else if(direction < .75 && zero-DIMENSION >= 0 && previous!=DIMENSION){
				move(-DIMENSION);
				previous = -DIMENSION;
				shifts++;
			}
			//If it is not on the bottom edge and doesn't undo
			else if (zero+DIMENSION < state.length && previous!=-DIMENSION){
				move(DIMENSION);
				previous = DIMENSION;
				shifts++;
			}
		}

		this.depth = 0; // Initialize the depth as 1

		//Calculates the two heuristic functions
		h1 = getHeuristic1(); 
		h2 = getHeuristic2();
		h3 = getHeuristic3();
	}


	/**
	 * Private helper method to generate initial game states by shifting the 0 tile.
	 * @param direction Direction to move the 0 tile.
	 */
	private void move(int direction){
		int temp = state[zero+direction];
		state[zero+direction] = 0;
		state[zero] = temp;
		zero = zero+direction;
		//System.out.println(this);
	}

	/**
	 * Private helper method to help find the optimal path to the goal state using the
	 * A* Search algorithm by adding children states to the priority queue.
	 * 
	 * @param pq The priority queue for the A* algorithm
	 */
	private int addChild(PriorityQueue<SlidingGameState> pq, HashSet<SlidingGameState> visited){

		int children = 0;
		if(zero%DIMENSION != 0){ // Empty is not on the left edge
			SlidingGameState left = child(-1);
			if(!visited.contains(left)){
				pq.add(left);
				children++;

			}
			//			children++;
		}

		if(zero%DIMENSION !=DIMENSION-1){ //Empty is not on the right edge
			SlidingGameState right = child(1);
			if(!visited.contains(right)){
				pq.add(right);
				children++;

			}
			//			children++;

		}

		if(zero-DIMENSION >= 0){ // Empty is not on the top edge
			SlidingGameState up = child(-DIMENSION);
			if(!visited.contains(up)){
				pq.add(up);
				children++;

			}
			//			children++;

		}

		if(zero+DIMENSION < state.length){ //Empty is not on the bottom edge
			SlidingGameState bot = child(DIMENSION);
			if(!visited.contains(bot)){
				pq.add(bot);
				children++;

			}
			//			children++;

		}
		return children;
	}

	/**
	 * Generates a child that results from a movement of the zero position 
	 * @param direction Direction of the movement of the empty zero position.
	 * @return A Game state with the zero position shifted.
	 */
	private SlidingGameState child(int direction){
		int[] child = new int[state.length];
		for(int i=0; i<state.length; i++){
			child[i] = state[i];
		}
		int temp = child[zero+direction];
		child[zero+direction] = 0;
		child[zero] = temp;
		return new SlidingGameState(child, this.depth+1, this.zero+direction);
	}

	/**
	 * Calculates the heuristic value of a tile as the number of tiles 
	 * misplaced from it's desired position.
	 * @return The number of positions is not correctly in the right 
	 * position on the board. 
	 */
	private int getHeuristic1(){
		int misplaced = 0;
		for (int i = 0; i < this.state.length; i++){
			if (this.state[i] != i){
				if(state[i]==0){
					zero = i;
					continue;
				}
				misplaced++;
			}
		}
		return misplaced;
	}

	/**
	 * Calculates the heuristic value as the number of moves for each tile to 
	 * reach their solved position if they could freely move.
	 * Can only be utilized if compareTo is changed from h1 to h2.
	 * @return The total number of moves to reach the goal state with free movement.
	 */
	private int getHeuristic2(){
		int misplaced = 0;
		for (int i = 0; i < this.state.length; i++){
			if (state[i] == 0){
				zero = i;
				continue;
			}
			misplaced += Math.abs(state[i] % DIMENSION - i % DIMENSION) + 
					Math.abs(state[i] / DIMENSION - i / DIMENSION);

		}
		return misplaced;
	}

	/**
	 * Calculates the heuristic value of Manhattan distance plus the number of additional moves required due
	 * to linear conflicts within rows and columns.
	 * 
	 * @return the Manhattan distance plus two times the number of linear conflicts.
	 */
	private int getHeuristic3(){
		int conf = 0;

		for(int i=0; i<DIMENSION; i++){
			int[] lc = new int[DIMENSION];
			int total = 0;

			//Rows
			for(int k=i*DIMENSION; k<i*DIMENSION+DIMENSION; k++){
				//If the current tile is on the right row and not zero
				if(state[k] / DIMENSION == i  && state[k]!=0){			
					for(int j = k+1; j<i*DIMENSION+DIMENSION; j++){
						//If the tiles after are on the right row, non-zero, and conflicting
						if(state[j] / DIMENSION == i  && state[j]<state[k] &&  state[j]!=0){
							//Count the linear conflicts in each
							lc[k-i*DIMENSION]++;
							lc[j-i*DIMENSION]++;
							total+=2;
						}
					}
				}
			}

			//Sum linear conflicts
			while(total>0){
				int largest = 0;
				for(int j=1; j<DIMENSION; j++){
					if(lc[j]>lc[largest]){
						largest = j;
					}
				}
				total-=lc[largest];
				lc[largest]=0;
				conf++;
				for(int j=0; j<DIMENSION; j++){
					if(lc[j]>0){
						lc[j]--;
						total--;
					}
				}
			}

			//Columns
			for(int k=i; k<DIMENSION*DIMENSION; k+=DIMENSION){
				//If the current tile is on the right column and not zero
				if(state[k] % DIMENSION == i  && state[k]!=0){
					for(int j = k+DIMENSION; j<DIMENSION*DIMENSION; j+=DIMENSION){
						//If the tiles after are on the right row, non-zero, and conflicting
						if(state[j] % DIMENSION == i  && state[j]<state[k] && state[j]!=0){
							lc[k/DIMENSION]++;
							lc[(j-k)/DIMENSION]++;
							total+=2;
						}
					}
				}
			}
			//Sum linear conflicts
			while(total>0){
				int largest = 0;
				for(int j=1; j<DIMENSION; j++){
					if(lc[j]>lc[largest]){
						largest = j;
					}
				}
				total-=lc[largest];
				lc[largest]=0;
				conf++;
				for(int j=0; j<DIMENSION; j++){
					if(lc[j]>0){
						lc[j]--;
						total--;
					}

				}
			}
		}
		return conf*2 + this.h2; //h2 is always calculated before h3 in our code
	}

	/**
	 * Compare a given SlidingGameBoard with current(this) object.
	 * compare their h values and total path costs.
	 * Allows PriorityQueues to correctly sort the GameStates based on 
	 * Path length and Heuristic values
	 */
	@SuppressWarnings("unused")
	@Override
	public int compareTo(SlidingGameState other) {
		if(HEURISTIC==1){
			int diff = (this.h1+this.depth) - (other.h1	 + other.depth);
			if(DEPTH && diff==0){
				return this.depth - other.depth;
			}
			if(HRESOLUTION && diff==0){
				return this.h1-other.h1;
			}
			return diff;
		}
		else if(HEURISTIC==2){
			int diff = (this.h2+this.depth) - (other.h2	 + other.depth);
			if(DEPTH && diff==0){
				return this.depth - other.depth;
			}
			if(HRESOLUTION && diff==0){
				return this.h2-other.h2;
			}
			return diff;
		}
		else{
			int diff = (this.h3+this.depth) - (other.h3	 + other.depth);
			if(DEPTH && diff==0){
				return this.depth - other.depth;

			}
			if(HRESOLUTION && diff==0){
				return this.h3-other.h3;
			}
			return diff;
		}
	}

	/**
	 * Override toString to to display the board as a square matrix
	 */
	@Override
	public String toString(){
		String s = "";
		for(int i=0; i<state.length; i+=DIMENSION){
			for(int j=i; j<i+DIMENSION; j++){
				s += state[j];
				s += "\t";
			}
			s+= "\n";
		}
		return s;
	}


	/**
	 * Overwrites the equals method so the contains method 
	 * recognizes repeated board states in the HashSet.
	 */
	@Override
	public boolean equals(Object obj){
		if(obj == this){
			return true;
		}
		if(obj == null || obj.getClass() != this.getClass()){
			return false;
		}
		SlidingGameState other = (SlidingGameState) obj;

		for(int i=0; i<this.state.length; i++){
			if(state[i]!=other.state[i]){
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns the hashCode for SlidingGameState object based on its state
	 * @return int - the hashcode is calculated by adding values of each state
	 * and multiply by a prime number
	 */
	@Override
	public int hashCode(){
		final int prime = 31;
		int hash = 1;
		for(int i=0; i<state.length; i++){
			hash *= prime;
			hash += state[i];
		}
		return hash;
	}

	/**
	 * Solves the game board from an initial state.
	 * @return An integer array with the depth of the solved state at index 0
	 * The ideal number of moves desired from the solved state at index 1
	 * The number of nodes expanded to reach the goal state at index 2.
	 */
	public int[] solveBoard(){ //Returns true if optimal
		PriorityQueue<SlidingGameState> pq = new PriorityQueue<SlidingGameState>();
		HashSet<SlidingGameState> visited = new HashSet<SlidingGameState>();
		pq.add(this);
		int steps = 1; //First node expanded
		while(!pq.isEmpty()){//Method to check if the board is solved
			SlidingGameState current = pq.poll();
			visited.add(current);
			//System.out.println("H1: "+current.h1+ " Depth: " + current.depth);
			//steps++;
			//			System.out.println(current);
			if(current.h1 == 0){
				//System.out.print(current);
				int[] results = new int[2];
				results[0] = current.depth;
				results[1] = steps;

				return results; 
			}
			steps += current.addChild(pq, visited);
		}
		return null; //Never unsolvable with our board generation method
	}
}
