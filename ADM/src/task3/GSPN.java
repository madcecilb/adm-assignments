package task3;

import java.util.ArrayList;

import task1.DTMC;
import task2.CTMC;

public class GSPN {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Double timeStep = 0.1;
		// TODO Auto-generated method stub
		ArrayList<ArrayList<Double>> trMat = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> state;
		for(int i = 0; i < 4; i++) {
			state = new ArrayList<Double>();
			if(i == 0) {
				state.add((-1.0/3.0) - (1.0/11.0));
				state.add(1.0/3.0);
				state.add(0.0);
				state.add(1.0/11.0);
			}
			if(i == 1) {
				state.add(0.0);
				state.add(-0.1 - (1.0/11.0));
				state.add(0.1);
				state.add(1.0/11.0);
			}
			if(i == 2) {
				state.add(0.0);
				state.add(0.0);
				state.add(-0.1);
				state.add(0.1);
			}
			if(i == 3) {
				state.add(0.0);
				state.add(0.0);
				state.add(0.0);
				state.add(0.0);
			}
			trMat.add(state);
		}
		
		trMat = CTMC.discretize(trMat, timeStep); 
		
		ArrayList<Double> initState = new ArrayList<Double>();
		initState.add(1.0); initState.add(0.0); initState.add(0.0); initState.add(0.0);
		
		ArrayList<ArrayList<Double>> states = DTMC.nSteps(trMat, initState, 12);
		//DTMC.printStates(DTMC.nSteps(trMat, initState, 0.0000001, 1000000));
		DTMC.printStates(states);
		
		System.out.println( 100 - 100*states.get(states.size() - 1).get(states.get(states.size() - 1).size() - 1));
		
	}

}
