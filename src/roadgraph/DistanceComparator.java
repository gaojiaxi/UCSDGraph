package roadgraph;

import java.util.Comparator;

public class DistanceComparator implements Comparator<MapNode>{
	public int compare(MapNode n1, MapNode n2) {
		if (n1.getDistanceToSource() < n2.getDistanceToSource()) {
			return -1;
		}
		else if (n1.getDistanceToSource() > n2.getDistanceToSource()) {
			return 1;
		}
		return 0;
	}
	
}
