package jminor;

import java.util.*;
import java.io.*;

public class MetroStationPlanner {

    private List<Map<String, Object>> data;
    private int[] clusters; // Store DBSCAN results
    private double[][] shortestPaths; // Store Floyd-Warshall results

    // Load dataset from a file
    public List<Map<String, Object>> loadData(String filePath) throws IOException {
        data = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String line;
        String[] headers = reader.readLine().split(",");
        
        while ((line = reader.readLine()) != null) {
            String[] values = line.split(",");
            Map<String, Object> row = new HashMap<>();
            for (int i = 0; i < headers.length; i++) {
                row.put(headers[i], values[i]);
            }
            data.add(row);
        }
        reader.close();
        return data;
    }

    // Getter for data
    public List<Map<String, Object>> getData() {
        return data;
    }

    // Preprocess data
    public List<Map<String, Object>> preprocessData(List<Map<String, Object>> data) {
        List<Map<String, Object>> filteredData = new ArrayList<>();
        for (Map<String, Object> row : data) {
            try {
                double latitude = Double.parseDouble(row.get("latitude").toString());
                double longitude = Double.parseDouble(row.get("longitude").toString());
                int population = Integer.parseInt(row.get("population").toString());
                if (latitude != 0 && longitude != 0 && population > 0) {
                    filteredData.add(row);
                }
            } catch (Exception e) {
                // Handle invalid data rows
            }
        }
        return filteredData;
    }

    // Calculate the geographical distance matrix
    public double[][] customGeographicalDistanceMatrix(List<Map<String, Object>> data) {
        int n = data.size();
        double[][] matrix = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double lat1 = Double.parseDouble(data.get(i).get("latitude").toString());
                    double lon1 = Double.parseDouble(data.get(i).get("longitude").toString());
                    double lat2 = Double.parseDouble(data.get(j).get("latitude").toString());
                    double lon2 = Double.parseDouble(data.get(j).get("longitude").toString());
                    matrix[i][j] = haversineDistance(lat1, lon1, lat2, lon2);
                } else {
                    matrix[i][j] = 0.0; // Distance to self is 0
                }
            }
        }
        return matrix;
    }

    // Haversine formula
    public double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
        final int R = 6371; // Radius of Earth in kilometers
        double latDistance = Math.toRadians(lat2 - lat1);
        double lonDistance = Math.toRadians(lon2 - lon1);
        double a = Math.sin(latDistance / 2) * Math.sin(latDistance / 2)
                + Math.cos(Math.toRadians(lat1)) * Math.cos(Math.toRadians(lat2))
                * Math.sin(lonDistance / 2) * Math.sin(lonDistance / 2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return R * c;
    }

    // Apply DBSCAN
    public int[] applyDBSCAN(double[][] distanceMatrix, double eps, int minPts) {
        int n = distanceMatrix.length;
        clusters = new int[n];
        Arrays.fill(clusters, -1); // Initialize all as noise (-1)

        int clusterId = 0;

        for (int i = 0; i < n; i++) {
            if (clusters[i] != -1) continue;
            List<Integer> neighbors = regionQuery(distanceMatrix, i, eps);
            if (neighbors.size() < minPts) {
                clusters[i] = -1;
            } else {
                expandCluster(distanceMatrix, clusters, i, neighbors, clusterId, eps, minPts);
                clusterId++;
            }
        }
        return clusters;
    }

    // Find neighbors
    private List<Integer> regionQuery(double[][] distanceMatrix, int point, double eps) {
        List<Integer> neighbors = new ArrayList<>();
        for (int i = 0; i < distanceMatrix.length; i++) {
            if (distanceMatrix[point][i] <= eps) {
                neighbors.add(i);
            }
        }
        return neighbors;
    }

    // Expand cluster
    private void expandCluster(double[][] distanceMatrix, int[] clusters, int point, List<Integer> neighbors,
                               int clusterId, double eps, int minPts) {
        clusters[point] = clusterId;
        Queue<Integer> queue = new LinkedList<>(neighbors);

        while (!queue.isEmpty()) {
            int current = queue.poll();
            if (clusters[current] == -1) {
                clusters[current] = clusterId;
            }
            if (clusters[current] != -1) continue;
            clusters[current] = clusterId;

            List<Integer> currentNeighbors = regionQuery(distanceMatrix, current, eps);
            if (currentNeighbors.size() >= minPts) {
                queue.addAll(currentNeighbors);
            }
        }
    }

    public static void main(String[] args) {
        System.out.println("Backend Ready!");
    }
}