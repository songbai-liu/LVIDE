package jmetal.metaheuristics.lsmoea;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import jmetal.core.SolutionSet;
import jmetal.metaheuristics.moead.Utils;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class GroupingStrategy {
	private int[][] group_; //Record the grouping results
	private int numOfGroups_; //number of groups
	private int numVar_; //number of variables 
	private int numObj_;//number of objectives
	
	public GroupingStrategy(int numOfGroup, int numVar, int numObj){
		numOfGroups_ = numOfGroup;
		numVar_ = numVar;
		group_ = new int[numOfGroups_][];
		numObj_ = numObj;
	}
	
	public int[][] linearGrouping(){
		int gSize_ = numVar_/numOfGroups_;
		for(int g=0;g<numOfGroups_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVar_-(numOfGroups_-1)*gSize_;//the variable size of the last group
		group_[numOfGroups_-1] = new int[lSize_];
		
		int t = 0;
		for(int g=0;g<numOfGroups_-1;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = t;
				t++;
			}
		}
		//assign variable to the last group
		for(int m=0;m<lSize_;m++){
			group_[numOfGroups_-1][m] = t;
			t++;
		}
		return group_;
	}
	
	public int[][] randomGrouping(){
		int gSize_ = numVar_/numOfGroups_;
		for(int g=0;g<numOfGroups_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVar_-(numOfGroups_-1)*gSize_;//the variable size of the last group
		group_[numOfGroups_-1] = new int[lSize_];
		int[] permutation = new int[numVar_];
		Utils.randomPermutation(permutation, numVar_);
		int t = 0;
		for(int g=0;g<numOfGroups_-1;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = permutation[t];
				t++;
			}
		}
		//assign variable to the last group
		for(int m=0;m<lSize_;m++){
			group_[numOfGroups_-1][m] = permutation[t];
			t++;
		}
		return group_;
	}
	
	public int[][] orderedGrouping(SolutionSet oldSet, SolutionSet newSet) throws JMException{
		int gSize_ = numVar_/numOfGroups_;
		for(int g=0;g<numOfGroups_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVar_-(numOfGroups_-1)*gSize_;//the variable size of the last group
		group_[numOfGroups_-1] = new int[lSize_];
		int solutionSize = oldSet.size();
		int[] index = new int[numVar_];
		double[] delta = new double[numVar_];
		for(int var=0;var<numVar_;var++){
			index[var] = var;
			delta[var] = 0.0;
		}
		for(int i=0;i<solutionSize;i++){
			XReal oldSol = new XReal(oldSet.get(i));
			XReal newSol = new XReal(newSet.get(i));
			for(int var=0;var<numVar_;var++){
				delta[var] += Math.abs(newSol.getValue(var)-oldSol.getValue(var))/(oldSol.getUpperBound(var)-oldSol.getLowerBound(var)); 
			}
		}
		for(int var=0;var<numVar_;var++){
			delta[var] = delta[var]/solutionSize;
		}
		QuickSort(delta, index, 0, numVar_-1);

		int t = 0;
		for(int g=0;g<numOfGroups_-1;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = index[t];
				t++;
			}
		}
		//assign variable to the last group
		for(int m=0;m<lSize_;m++){
			group_[numOfGroups_-1][m] = index[t];
			t++;
		}
		return group_;
	}
	
	public int[][] contributionGrouping(double[][] contribution){
		int gSize_ = numVar_/numOfGroups_;
		for(int g=0;g<numOfGroups_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVar_-(numOfGroups_-1)*gSize_;//the variable size of the last group
		group_[numOfGroups_-1] = new int[lSize_];
		// dominateMe[i] contains the number of solutions dominating i
		int[] dominateMe = new int[numVar_];
		// iDominate[k] contains the list of solutions dominated by k
		List<Integer>[] iDominate = new List[numVar_];
		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[numVar_+1];
		// flagDominate is an auxiliar encodings.variable
		int flagDominate;
		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();
		for (int d = 0; d < numVar_; d++) {
			// Initialize the list of individuals that i dominate and the number
			// of individuals that dominate me
			iDominate[d] = new LinkedList<Integer>();
			dominateMe[d] = 0;
		}			
		for (int p = 0; p < (numVar_ - 1); p++) {
			for (int q = p + 1; q < numVar_; q++) {
				flagDominate = dominance_Compare(contribution[p],contribution[q]);
				if (flagDominate == -1) {
					iDominate[p].add(q);
					dominateMe[q]++;
				} else if (flagDominate == 1) {
					iDominate[q].add(p);
					dominateMe[p]++;
				}
			}
		}
		for (int p = 0; p < numVar_; p++) {
			if (dominateMe[p] == 0) {
				front[0].add(p);
			}
		}
		// Obtain the rest of fronts
		int i = 0;
		Iterator<Integer> it1, it2; // Iterators
		while (front[i].size() != 0) {
			i++;
			it1 = front[i - 1].iterator();
			while (it1.hasNext()) {
				it2 = iDominate[it1.next()].iterator();
				while (it2.hasNext()) {
					int index = it2.next();
					dominateMe[index]--;
					if (dominateMe[index] == 0) {
						front[i].add(index);
					}
				}
			}
		}
		//Contribution-based grouping
		int L = front.length;
		int[] tt = new int[numVar_];
		int u = 0;
		for(int l=0;l<L;l++){
			for(int t=0;t<front[l].size();t++){
				tt[u] = front[l].get(t);
				u++;
			}
		}
		int t = 0;
		for(int g=0;g<numOfGroups_-1;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = tt[t];
				t++;
			}
		}
		//assign variable to the last group
		for(int m=0;m<lSize_;m++){
			group_[numOfGroups_-1][m] = tt[t];
			t++;
		}		
		return group_;
	}
	
	public int[][] contributionGrouping1(double[][] contribution){
		//int gSize_ = numVar_/numOfGroups_;
		/*for(int g=0;g<numOfGroups_-1;g++){
			group_[g] = new int[gSize_];
		}
		int lSize_ = numVar_-(numOfGroups_-1)*gSize_;//the variable size of the last group
		group_[numOfGroups_-1] = new int[lSize_];*/
		// dominateMe[i] contains the number of solutions dominating i
		int[] dominateMe = new int[numVar_];
		// iDominate[k] contains the list of solutions dominated by k
		List<Integer>[] iDominate = new List[numVar_];
		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[numVar_+1];
		// flagDominate is an auxiliar encodings.variable
		int flagDominate;
		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();
		for (int d = 0; d < numVar_; d++) {
			// Initialize the list of individuals that i dominate and the number
			// of individuals that dominate me
			iDominate[d] = new LinkedList<Integer>();
			dominateMe[d] = 0;
		}			
		for (int p = 0; p < (numVar_ - 1); p++) {
			for (int q = p + 1; q < numVar_; q++) {
				flagDominate = dominance_Compare(contribution[p],contribution[q]);
				if (flagDominate == -1) {
					iDominate[p].add(q);
					dominateMe[q]++;
				} else if (flagDominate == 1) {
					iDominate[q].add(p);
					dominateMe[p]++;
				}
			}
		}
		for (int p = 0; p < numVar_; p++) {
			if (dominateMe[p] == 0) {
				front[0].add(p);
			}
		}
		// Obtain the rest of fronts
		int i = 0;
		Iterator<Integer> it1, it2; // Iterators
		while (front[i].size() != 0) {
			i++;
			it1 = front[i - 1].iterator();
			while (it1.hasNext()) {
				it2 = iDominate[it1.next()].iterator();
				while (it2.hasNext()) {
					int index = it2.next();
					dominateMe[index]--;
					if (dominateMe[index] == 0) {
						front[i].add(index);
					}
				}
			}
		}
		//Contribution-based grouping
		int L = front.length;
		int[][] group = new int[L][];
		for(int l=0;l<L;l++){
			group[l] = new int[front[l].size()];
			int u = 0;
			for(int t=0;t<front[l].size();t++){
				group[l][u] = front[l].get(t);
				u++;
			}
		}
		return group;
	}
	
	public int[][] broadeningGrouping(int[][] group){
		int numOfGroup = group.length;
		int[][] groups = new int[numOfGroup][];
		int T_ = 0;
		groups[0] = group[0];
		while(T_ < numOfGroup-1){
			T_++;
			int[] groupA = groups[T_-1];
			int[] groupB = group[T_];
			groups[T_] = new int[groupA.length + groupB.length];
			System.arraycopy(groupA, 0, groups[T_], 0, groupA.length);
			System.arraycopy(groupB, 0, groups[T_], groupA.length, groupB.length);
		}
		return groups;
	}
	
	static void QuickSort(double[] array, int[] idx, int from, int to) {
		if (from < to) {
			double temp = array[to];
			int tempIdx = idx[to];
			int i = from - 1;
			for (int j = from; j < to; j++) {
				if (array[j] <= temp) {
					i++;
					double tempValue = array[j];
					array[j] = array[i];
					array[i] = tempValue;
					int tempIndex = idx[j];
					idx[j] = idx[i];
					idx[i] = tempIndex;
				}
			}
			array[to] = array[i + 1];
			array[i + 1] = temp;
			idx[to] = idx[i + 1];
			idx[i + 1] = tempIdx;
			QuickSort(array, idx, from, i);
			QuickSort(array, idx, i + 1, to);
		}
	}
	
	public int dominance_Compare(double[] value1, double[] value2){
		//int numObj_ = problem_.getNumberOfObjectives();
		if (value1==null)
		    return 1;
		else if (value2 == null)
		    return -1;

		int dominate1 ; // dominate1 indicates if some objective of solution1 
		                    // dominates the same objective in solution2. dominate2
		int dominate2 ; // is the complementary of dominate1.

		dominate1 = 0 ; 
		dominate2 = 0 ;
		    
		int flag; //stores the result of the comparison
		    
		// Equal number of violated constraints. Applying a dominance Test then
		double fit1, fit2;
		for (int i = 0; i < numObj_; i++) {
		   fit1 = value1[i];
		   fit2 = value2[i];
		   if (fit1 > fit2) {
		      flag = -1;
		   } else if (fit1 < fit2) {
		      flag = 1;
		   } else {
		      flag = 0;
		   }
		      
		   if (flag == -1) {
		      dominate1 = 1;
		   }
		      
		   if (flag == 1) {
		      dominate2 = 1;           
		   }
		}
		            
		if (dominate1 == dominate2) {            
		   return 0; //No one dominate the other
		}
		if (dominate1 == 1) {
		   return -1; // solution1 dominate
		}
		return 1;    // solution2 dominate  
	}//dominance_Compare

}
