from kneed import KneeLocator
import matplotlib.pyplot as plt
import json

def plot_pareto_front(x, y, s, curve, direction, x_label):
    knee = KneeLocator(x, y, s, curve=curve, direction=direction)
    (knee_x, knee_y) = (knee.knee, knee.knee_y)
    
    print(knee_x, knee_y)
    plt.plot(x, y, marker='o', linewidth=2, markersize=1)
    plt.ylabel(x_label)
    plt.plot([knee_x], [knee_y], marker='*', color ='red', markersize=10)
    plt.xlabel("Regularization Weight")
    plt.title(f"Regularized Parsimony Pareto Front (s={s})")
    plt.show()

def parse_results(dir, weights):
    parsimony_score = []
    num_edges = []
    obj = []
    for weight in weights:
        fname = f"{dir}/inferred_labeling_regularized_w={weight}_results.json"
        with open(fname, "r") as file:
            data = json.load(file)
            obj.append(data["objective"])
            parsimony_score.append(data["migrations"])
            num_edges.append(data["migration_pattern_edges"])

    return obj, parsimony_score, num_edges

if __name__ == "__main__":
    dir = "./examples/cancer_evolution/simulations/calibrate_regularization_small"
    weights = [0, 0.001, 0.01, 0.1, 0.3, 0.5, 0.8, 1, 2, 3]
    (obj, parsimony_score, num_edges) = parse_results(dir, weights)

    print(parsimony_score)
    print(num_edges)

    plot_pareto_front(weights, parsimony_score, 0.5, "convex", "increasing", "Parsimony")
    plot_pareto_front(weights, parsimony_score, 1, "convex", "increasing", "Parsimony")
    plot_pareto_front(weights, num_edges, 1, "concave", "decreasing", "Edges in Migration pattern")

    print("This file will run multiple iterations of ancestral reconstruction defining a pareto front then find the knee/elbow")
    exit