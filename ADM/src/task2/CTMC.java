

package task2;

import java.util.ArrayList;

import task1.DTMC;

public class CTMC {
	
	public static ArrayList<ArrayList<Double>> discretize(ArrayList<ArrayList<Double>> oldMatrix, Double timeStep) {
		ArrayList<ArrayList<Double>> newMatrix = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> row;
		for(int i = 0; i < oldMatrix.size(); i++) {
			row = new ArrayList<Double>();
			for(int j = 0; j < oldMatrix.get(i).size(); j++) {
				if( i == j ) {
					row.add(1 + timeStep*oldMatrix.get(i).get(j));
				}
				else{
					row.add(timeStep*oldMatrix.get(i).get(j));
				}
				
			}
			newMatrix.add(row);
		}
		
		return newMatrix;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Double timeStep = 0.1;
		// TODO Auto-generated method stub
		ArrayList<ArrayList<Double>> trMat = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> state;
		for(int i = 0; i < 4; i++) {
			state = new ArrayList<Double>();
			if(i == 0) {
				state.add(-0.7);
				state.add(0.6);
				state.add(0.0);
				state.add(0.1);
			}
			if(i == 1) {
				state.add(0.0);
				state.add(-0.3);
				state.add(0.2);
				state.add(0.1);
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
		
		trMat = discretize(trMat, timeStep); 
		
		ArrayList<Double> initState = new ArrayList<Double>();
		initState.add(1.0); initState.add(0.0); initState.add(0.0); initState.add(0.0);
		
		DTMC.printStates(DTMC.nSteps(trMat, initState, 0.000000001, 1000000));
	}

}
