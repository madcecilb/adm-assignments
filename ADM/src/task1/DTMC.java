package task1;

import java.lang.reflect.Array;
import java.util.ArrayList;

public class DTMC {
	
	private static Double stateDifference(ArrayList<Double> curStep, ArrayList<Double> prevStep) {
		Double diff = 0.0;
		for(int i = 0; i < prevStep.size(); i++) {
			diff += Math.pow(prevStep.get(i) - curStep.get(i), 2);
		}
		return diff;
	}
	
	public static ArrayList<Double> nextStep(ArrayList<ArrayList<Double>> transition_matrix, ArrayList<Double> state) {
		ArrayList<Double> newState = new ArrayList<Double>();
		for(int i = 0; i < transition_matrix.size(); i++) {
			for(int j = 0; j < state.size(); j++ ) {
				if( i == 0) {
					newState.add(state.get(i) * transition_matrix.get(i).get(j));
				}
				else{
					newState.set(j, newState.get(j) + state.get(i) * transition_matrix.get(i).get(j));
				}
			}
		}

		return newState;		
	}

	public static ArrayList<ArrayList<Double>> nSteps(ArrayList<ArrayList<Double>> transition_matrix, 
			ArrayList<Double> initialStep, int numberOfSteps) {
		ArrayList<ArrayList<Double>> states = new ArrayList<ArrayList<Double>>();
		states.add(initialStep);
		ArrayList<Double> currentStep = initialStep;
		for(int i = 0; i < numberOfSteps; i++) {
			currentStep = nextStep(transition_matrix, currentStep);
			states.add(currentStep);
		}
		
		return states;		
	}
	
	public static ArrayList<ArrayList<Double>> nSteps(ArrayList<ArrayList<Double>> transition_matrix, 
			ArrayList<Double> initialStep, double threshold, int stepLimit) {
		ArrayList<ArrayList<Double>> states = new ArrayList<ArrayList<Double>>();
		states.add(initialStep);
		ArrayList<Double> currentStep = initialStep;
		ArrayList<Double> prevStep;
		for(int i = 0; i < stepLimit; i++) {
			prevStep = currentStep;
			currentStep = nextStep(transition_matrix, prevStep);
			states.add(currentStep);
			if(stateDifference(currentStep, prevStep) < threshold) {
				break;
			}			
		}
		
		return states;		
	}
	
	public static void printStates(ArrayList<ArrayList<Double>> steps) {
		for(int i = 1; i <= steps.size(); i++ ){
			System.out.println(i + "------------------------//" );
			for (Double val : steps.get(i -1)) {
				System.out.print(val + " ; ");
			}
			System.out.println();
		}
	}

	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ArrayList<ArrayList<Double>> trMat = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> state;
		for(int i = 0; i < 4; i++) {
			state = new ArrayList<Double>();
			if(i == 0) {
				state.add(0.6);
				state.add(0.3);
				state.add(0.0);
				state.add(0.1);
			}
			if(i == 1) {
				state.add(0.0);
				state.add(0.75);
				state.add(0.15);
				state.add(0.1);
			}
			if(i == 2) {
				state.add(0.0);
				state.add(0.0);
				state.add(0.9);
				state.add(0.1);
			}
			if(i == 3) {
				state.add(0.0);
				state.add(0.0);
				state.add(0.0);
				state.add(1.0);
			}
			trMat.add(state);
		}
		
		ArrayList<Double> initState = new ArrayList<Double>();
		initState.add(1.0); initState.add(0.0); initState.add(0.0); initState.add(0.0);
		
		printStates(nSteps(trMat, initState, 0.000001, 1000));
	}

}
