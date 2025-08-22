import React, { useState, useEffect, useRef, useCallback } from 'react';

const AlgorithmVisualizer = () => {
  // State variables
  const [selectedAlgorithm, setSelectedAlgorithm] = useState('sort:bubble');
  const [speed, setSpeed] = useState(50);
  const [count, setCount] = useState(5);
  const [target, setTarget] = useState('');
  const [customArray, setCustomArray] = useState('');
  const [activeTab, setActiveTab] = useState('sorting');
  const [array, setArray] = useState([]);
  const [arrayElements, setArrayElements] = useState([]);
  const [graph, setGraph] = useState({});
  const [graphNodes, setGraphNodes] = useState([]);
  const [graphEdges, setGraphEdges] = useState([]);
  const [startNode, setStartNode] = useState('');
  const [traversalSteps, setTraversalSteps] = useState([]);
  const [currentStep, setCurrentStep] = useState(0);
  const [comparisons, setComparisons] = useState(0);
  const [swaps, setSwaps] = useState(0);
  const [iterations, setIterations] = useState(0);
  const [isRunning, setIsRunning] = useState(false);
  const [searchResult, setSearchResult] = useState(null);
  const [isPaused, setIsPaused] = useState(false);
  const [graphInput, setGraphInput] = useState('');
  const [graphType, setGraphType] = useState('directed');
  const [edgeWeight, setEdgeWeight] = useState(1);
  const [fromNode, setFromNode] = useState('');
  const [toNode, setToNode] = useState('');
  const [nodePositions, setNodePositions] = useState({});
  const [nodeValue, setNodeValue] = useState('');

  // Refs
  const abortControllerRef = useRef(null);
  const autoPlayIntervalRef = useRef(null);
  const arrayContainerRef = useRef(null);
  const stepsContainerRef = useRef(null);
  const graphContainerRef = useRef(null);
  const stepRefs = useRef([]);

  // Algorithm code examples
  const algorithmCode = {
    'bubble': `function bubbleSort(arr) { 
  let n = arr.length; 
  for (let i = 0; i < n-1; i++) { 
    for (let j = 0; j < n-i-1; j++) { 
      if (arr[j] > arr[j+1]) { 
        // Swap arr[j] and arr[j+1] 
        let temp = arr[j]; 
        arr[j] = arr[j+1]; 
        arr[j+1] = temp; 
      } 
    } 
  } 
}`,

    'selection': `function selectionSort(arr) { 
  let n = arr.length; 
  for (let i = 0; i < n-1; i++) { 
    let minIdx = i; 
    for (let j = i+1; j < n; j++) { 
      if (arr[j] < arr[minIdx]) { 
        minIdx = j; 
      } 
    } 
    let temp = arr[minIdx]; 
    arr[minIdx] = arr[i]; 
    arr[i] = temp; 
  } 
}`,

    'insertion': `function insertionSort(arr) { 
  let n = arr.length; 
  for (let i = 1; i < n; i++) { 
    let key = arr[i]; 
    let j = i - 1; 
    while (j >= 0 && arr[j] > key) { 
      arr[j+1] = arr[j]; 
      j = j - 1; 
    } 
    arr[j+1] = key; 
  } 
}`,

    'merge': `function mergeSort(arr) {
  if (arr.length <= 1) return arr;

  let mid = Math.floor(arr.length / 2);
  let left = mergeSort(arr.slice(0, mid));
  let right = mergeSort(arr.slice(mid));

  return merge(left, right);
}

function merge(left, right) {
  let result = [];
  let i = 0, j = 0;

  while (i < left.length && j < right.length) {
    if (left[i] < right[j]) {
      result.push(left[i]);
      i++;
    } else {
      result.push(right[j]);
      j++;
    }
  }

  return result.concat(left.slice(i)).concat(right.slice(j));
}`,

    'quick': `function quickSort(arr, low = 0, high = arr.length - 1) {
  if (low < high) {
    let pi = partition(arr, low, high);
    quickSort(arr, low, pi - 1);
    quickSort(arr, pi + 1, high);
  }
}

function partition(arr, low, high) {
  let pivot = arr[high];
  let i = low - 1;

  for (let j = low; j < high; j++) {
    if (arr[j] < pivot) {
      i++;
      [arr[i], arr[j]] = [arr[j], arr[i]];
    }
  }

  [arr[i + 1], arr[high]] = [arr[high], arr[i + 1]];
  return i + 1;
}`,

    'heap': `function heapSort(arr) {
  let n = arr.length;

  // Build max heap
  for (let i = Math.floor(n / 2) - 1; i >= 0; i--) {
    heapify(arr, n, i);
  }

  // Extract elements from heap
  for (let i = n - 1; i > 0; i--) {
    [arr[0], arr[i]] = [arr[i], arr[0]];
    heapify(arr, i, 0);
  }
}

function heapify(arr, n, i) {
  let largest = i;
  let left = 2 * i + 1;
  let right = 2 * i + 2;

  if (left < n && arr[left] > arr[largest]) {
    largest = left;
  }

  if (right < n && arr[right] > arr[largest]) {
    largest = right;
  }

  if (largest !== i) {
    [arr[i], arr[largest]] = [arr[largest], arr[i]];
    heapify(arr, n, largest);
  }
}`,

    'linear': `function linearSearch(arr, target) {
  for (let i = 0; i < arr.length; i++) {
    if (arr[i] === target) {
      return i; // Found at index i
    }
  }
  return -1; // Not found
}`,

    'binary': `function binarySearch(arr, target) {
  let left = 0;
  let right = arr.length - 1;
  while (left <= right) {
    let mid = Math.floor((left + right) / 2);
    if (arr[mid] === target) {
      return mid;
    } else if (arr[mid] < target) {
      left = mid + 1;
    } else {
      right = mid - 1;
    }
  }
  return -1;
}`,

    'jump': `function jumpSearch(arr, target) {
  let n = arr.length;
  let step = Math.floor(Math.sqrt(n));
  let prev = 0;

  while (arr[Math.min(step, n) - 1] < target) {
    prev = step;
    step += Math.floor(Math.sqrt(n));
    if (prev >= n) return -1;
  }

  while (arr[prev] < target) {
    prev++;
    if (prev == Math.min(step, n)) return -1;
  }

  if (arr[prev] == target) return prev;
  return -1;
}`,

    'interpolation': `function interpolationSearch(arr, target) {
  let low = 0;
  let high = arr.length - 1;

  while (low <= high && target >= arr[low] && target <= arr[high]) {
    if (low === high) {
      if (arr[low] === target) return low;
      return -1;
    }

    let pos = low + Math.floor(((target - arr[low]) * (high - low)) / (arr[high] - arr[low]));

    if (arr[pos] === target) return pos;
    if (arr[pos] < target) low = pos + 1;
    else high = pos - 1;
  }
  return -1;
}`,

    'dijkstra': `function dijkstra(graph, start) {
  let distances = {};
  let visited = new Set();

  // Initialize distances
  for (let vertex in graph) {
    distances[vertex] = Infinity;
  }
  distances[start] = 0;

  while (true) {
    let minVertex = null;
    for (let vertex in graph) {
      if (!visited.has(vertex) &&
          (minVertex === null || distances[vertex] < distances[minVertex])) {
        minVertex = vertex;
      }
    }

    if (minVertex === null) break;
    visited.add(minVertex);

    for (let neighbor in graph[minVertex]) {
      let distance = distances[minVertex] + graph[minVertex][neighbor];
      if (distance < distances[neighbor]) {
        distances[neighbor] = distance;
      }
    }
  }
  return distances;
}`,

    'bellmanford': `function bellmanFord(graph, start) {
  let distances = {};
  for (let vertex in graph) {
    distances[vertex] = Infinity;
  }
  distances[start] = 0;

  for (let i = 0; i < Object.keys(graph).length - 1; i++) {
    for (let u in graph) {
      for (let v in graph[u]) {
        if (distances[u] + graph[u][v] < distances[v]) {
          distances[v] = distances[u] + graph[u][v];
        }
      }
    }
  }

  // Check for negative weight cycles
  for (let u in graph) {
    for (let v in graph[u]) {
      if (distances[u] + graph[u][v] < distances[v]) {
        return "Graph contains a negative weight cycle";
      }
    }
  }

  return distances;
}`,

    'prims': `function prims(graph) {
  let result = [];
  let visited = new Set();
  if (Object.keys(graph).length === 0) return result;

  let startNode = Object.keys(graph)[0];
  visited.add(startNode);

  while (visited.size < Object.keys(graph).length) {
    let minEdge = null;
    for (let u of visited) {
      for (let v in graph[u]) {
        if (!visited.has(v)) {
          if (minEdge === null || graph[u][v] < minEdge.weight) {
            minEdge = { u, v, weight: graph[u][v] };
          }
        }
      }
    }
    if (minEdge) {
      result.push(minEdge);
      visited.add(minEdge.v);
    } else {
      break; // No more edges to add
    }
  }
  return result;
}`,

    'kruskal': `function kruskal(graph) {
  let edges = [];
  for (let u in graph) {
    for (let v in graph[u]) {
      edges.push({ u, v, weight: graph[u][v] });
    }
  }
  edges.sort((a, b) => a.weight - b.weight);

  let parent = {};
  Object.keys(graph).forEach(node => parent[node] = node);

  function find(i) {
    if (parent[i] === i) return i;
    return find(parent[i]);
  }

  function union(i, j) {
    let rootI = find(i);
    let rootJ = find(j);
    if (rootI !== rootJ) {
      parent[rootJ] = rootI;
      return true;
    }
    return false;
  }

  let mst = [];
  for (let edge of edges) {
    if (union(edge.u, edge.v)) {
      mst.push(edge);
    }
  }
  return mst;
}`,

    'graph:bfs': `function graphBFS(graph, start) {
  let visited = new Set();
  let queue = [start];
  let result = [];

  visited.add(start);

  while (queue.length > 0) {
    let vertex = queue.shift();
    result.push(vertex);

    for (let neighbor of graph[vertex]) {
      if (!visited.has(neighbor)) {
        visited.add(neighbor);
        queue.push(neighbor);
      }
    }
  }
  return result;
}`,

    'graph:dfs': `function graphDFS(graph, start) {
  let visited = new Set();
  let result = [];

  function traverse(vertex) {
    visited.add(vertex);
    result.push(vertex);

    for (let neighbor of graph[vertex]) {
      if (!visited.has(neighbor)) {
        traverse(neighbor);
      }
    }
  }
  traverse(start);
  return result;
}`
  };

  // Initialize with random array
  useEffect(() => {
    generateArray();
    initializeGraph();
  }, []);

  // Update array elements when array changes
  useEffect(() => {
    const elements = array.map((value, index) => ({
      value,
      index,
      state: 'default',
      isTarget: false
    }));
    setArrayElements(elements);
  }, [array]);

  // Auto-scroll to current step - improved version
  useEffect(() => {
    if (stepRefs.current[currentStep]) {
      stepRefs.current[currentStep].scrollIntoView({
        behavior: 'smooth',
        block: 'center',
      });
    }
  }, [currentStep, traversalSteps.length]);

  // Initialize graph with some nodes
  const initializeGraph = () => {
    const initialGraph = {
      'A': { 'B': 4, 'C': 2 },
      'B': { 'A': 4, 'C': 1, 'D': 5 },
      'C': { 'A': 2, 'B': 1, 'D': 8, 'E': 10 },
      'D': { 'B': 5, 'C': 8, 'E': 2, 'F': 6 },
      'E': { 'C': 10, 'D': 2, 'F': 3 },
      'F': { 'D': 6, 'E': 3 }
    };
    setGraph(initialGraph);
    const nodes = Object.keys(initialGraph);
    setGraphNodes(nodes);

    const edges = [];
    for (let node in initialGraph) {
      for (let neighbor in initialGraph[node]) {
        edges.push({
          from: node,
          to: neighbor,
          weight: initialGraph[node][neighbor]
        });
      }
    }
    setGraphEdges(edges);
    setStartNode('A');
  };

  // Generate random array
  const generateArray = () => {
    const newArray = Array.from({ length: count }, () => Math.floor(Math.random() * 99) + 1);
    setArray(newArray);
    setSearchResult(null);
    resetStats();
  };

  // Apply custom array
  const applyCustomArray = () => {
    if (!customArray) return;
    const newArray = customArray.split(',')
      .map(val => parseInt(val.trim()))
      .filter(val => !isNaN(val));
      
    if (newArray.length === 0) {
      alert('Please enter valid numbers separated by commas');
      return;
    }

    if (newArray.length > 20) {
      alert('Array truncated to 20 elements for better visualization');
      newArray.length = 20;
    }

    setArray(newArray);
    setSearchResult(null);
    resetStats();
  };

  // Parse graph input
  const parseGraphInput = () => {
    if (!graphInput) return;
    try {
      const parsedGraph = JSON.parse(graphInput);
      setGraph(parsedGraph);
      setGraphNodes(Object.keys(parsedGraph));
      
      const edges = [];
      for (let node in parsedGraph) {
        for (let neighbor in parsedGraph[node]) {
          edges.push({
            from: node,
            to: neighbor,
            weight: parsedGraph[node][neighbor]
          });
        }
      }
      setGraphEdges(edges);
      
      if (!startNode || !Object.keys(parsedGraph).includes(startNode)) {
        setStartNode(Object.keys(parsedGraph)[0] || '');
      }
    } catch (error) {
      alert('Invalid graph format. Please use valid JSON format: {"A": {"B": 1}, "B": {"A": 1}}');
    }
  };

  // Add a node to the graph
  const addNode = () => {
    if (!nodeValue || graphNodes.includes(nodeValue)) return;
    const newGraph = { ...graph };
    newGraph[nodeValue] = {};
    setGraph(newGraph);
    setGraphNodes([...graphNodes, nodeValue]);
    setNodeValue('');
  };

  // Add an edge to the graph
  const addEdge = () => {
    if (!fromNode || !toNode || fromNode === toNode) return;
    const newGraph = { ...graph };
    if (!newGraph[fromNode]) newGraph[fromNode] = {};
    if (!newGraph[toNode]) newGraph[toNode] = {};

    newGraph[fromNode][toNode] = edgeWeight;
    if (graphType === 'undirected') {
      newGraph[toNode][fromNode] = edgeWeight;
    }

    setGraph(newGraph);

    // Update edges
    const edges = [];
    for (let node in newGraph) {
      for (let neighbor in newGraph[node]) {
        edges.push({
          from: node,
          to: neighbor,
          weight: newGraph[node][neighbor]
        });
      }
    }
    setGraphEdges(edges);

    setFromNode('');
    setToNode('');
    setEdgeWeight(1);
  };

  // Remove a node from the graph
  const removeNode = (node) => {
    const newGraph = { ...graph };
    delete newGraph[node];
    // Remove all edges to this node
    for (let n in newGraph) {
      if (newGraph[n][node]) {
        delete newGraph[n][node];
      }
    }

    setGraph(newGraph);
    setGraphNodes(graphNodes.filter(n => n !== node));

    if (startNode === node) {
      setStartNode(graphNodes.length > 0 ? graphNodes[0] : '');
    }
  };

  // Remove an edge from the graph
  const removeEdge = (from, to) => {
    const newGraph = { ...graph };
    if (newGraph[from] && newGraph[from][to]) {
      delete newGraph[from][to];
    }
    if (graphType === 'undirected' && newGraph[to] && newGraph[to][from]) {
      delete newGraph[to][from];
    }
    setGraph(newGraph);

    // Update edges
    const edges = [];
    for (let node in newGraph) {
      for (let neighbor in newGraph[node]) {
        edges.push({
          from: node,
          to: neighbor,
          weight: newGraph[node][neighbor]
        });
      }
    }
    setGraphEdges(edges);
  };

  // Reset statistics
  const resetStats = () => {
    setComparisons(0);
    setSwaps(0);
    setIterations(0);
    setTraversalSteps([]);
    setCurrentStep(0);
  };

  // Sleep function for animation delays
  const sleep = () => {
    return new Promise(resolve => setTimeout(resolve, 210 - speed));
  };

  // Add a traversal step with enhanced highlighting
  const addTraversalStep = (stepInfo) => {
    const { description, highlightIndices = [], stateChanges = {}, arrayState, graphState, comparisons, swaps, iterations } = stepInfo;
    const newStep = {
      description,
      highlightIndices,
      stateChanges,
      comparisons,
      swaps,
      iterations,
      arrayState: arrayState ? [...arrayState] : null,
      graphState: graphState ? JSON.parse(JSON.stringify(graphState)) : null,
    };

    setTraversalSteps(prev => {
      const newSteps = [...prev, newStep];
      // Auto-advance current step when running
      if (isRunning) {
        setCurrentStep(newSteps.length - 1);
      }
      return newSteps;
    });
  };

  // Update array element states with enhanced visual feedback
  const updateElementStates = (indices, state) => {
    setArrayElements(prev => prev.map((el, i) => {
      if (indices.includes(i)) {
        return {...el, state};
      }
      return el;
    }));
  };

  // Update specific element state
  const updateElementState = (index, state) => {
    setArrayElements(prev => prev.map((el, i) =>
      i === index ? {...el, state} : el
    ));
  };

  // Reset all element states
  const resetElementStates = () => {
    setArrayElements(prev => prev.map(el => ({...el, state: 'default'})));
  };

  // Start algorithm visualization
  const runAlgorithm = async () => {
    if (isRunning) return;
    setIsRunning(true);
    resetStats();
    resetElementStates();

    // Create a new AbortController for this run
    abortControllerRef.current = new AbortController();
    const signal = abortControllerRef.current.signal;

    try {
      const [group, name] = selectedAlgorithm.split(':');
      
      if (group === 'sort') {
        await runSortingAlgorithm(name, signal);
      } else if (group === 'search') {
        await runSearchAlgorithm(name, signal);
      } else if (group === 'graph') {
        await runGraphAlgorithm(name, signal);
      }
    } catch (error) {
      if (error.name !== 'AbortError') {
        console.error('Algorithm execution error:', error);
      }
    } finally {
      setIsRunning(false);
    }
  };

  // Stop algorithm execution
  const stopAlgorithm = () => {
    if (abortControllerRef.current) {
      abortControllerRef.current.abort();
    }
    if (autoPlayIntervalRef.current) {
      clearInterval(autoPlayIntervalRef.current);
      autoPlayIntervalRef.current = null;
    }

    setIsRunning(false);
    setIsPaused(false);
  };

  // Pause/Resume algorithm execution
  const togglePause = () => {
    setIsPaused(!isPaused);
  };

  // Enhanced Bubble Sort implementation with better highlighting
  const bubbleSort = async (signal) => {
    let arr = [...array];
    let n = arr.length;
    let localComparisons = 0, localSwaps = 0, localIterations = 0;
    addTraversalStep({ description: "🚀 Starting Bubble Sort - comparing adjacent elements", arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });

    for (let i = 0; i < n - 1; i++) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(s => s + 1);
      addTraversalStep({ description: `🔄 Pass ${i+1}/${n-1} - largest element will bubble to position ${n-i-1}`, arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      
      for (let j = 0; j < n - i - 1; j++) {
        if (signal.aborted) return;
        while (isPaused) await sleep();
        
        updateElementStates([j, j+1], 'active');
        localComparisons++;
        setComparisons(c => c + 1);
        addTraversalStep({ description: `🔍 Comparing ${arr[j]} and ${arr[j+1]} at positions ${j} and ${j+1}`, highlightIndices: [j, j+1], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
        
        await sleep();
        
        if (arr[j] > arr[j+1]) {
          updateElementStates([j, j+1], 'swap');
          addTraversalStep({ description: `🔄 Swapping ${arr[j]} and ${arr[j+1]} (${arr[j]} > ${arr[j+1]})`, highlightIndices: [j, j+1], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
          
          [arr[j], arr[j+1]] = [arr[j+1], arr[j]];
          setArray([...arr]);
          localSwaps++;
          setSwaps(s => s + 1);
          
          await sleep();
        } else {
          addTraversalStep({ description: `✅ No swap needed (${arr[j]} ≤ ${arr[j+1]})`, highlightIndices: [j, j+1], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
        }
        
        updateElementStates([j, j+1], 'default');
      }
      
      updateElementState(n - i - 1, 'ok');
      addTraversalStep({ description: `🎯 Element ${arr[n-i-1]} is now in its final position (index ${n-i-1})`, highlightIndices: [n-i-1], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
    }

    updateElementState(0, 'ok');
    addTraversalStep({ description: "🎉 Bubble Sort completed! All elements are now sorted.", arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
  };

  // Enhanced Selection Sort implementation
  const selectionSort = async (signal) => {
    let arr = [...array];
    let n = arr.length;
    let localComparisons = 0, localSwaps = 0, localIterations = 0;
    addTraversalStep({ description: "🚀 Starting Selection Sort - finding minimum elements", arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });

    for (let i = 0; i < n - 1; i++) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(s => s + 1);
      let minIdx = i;
      
      updateElementState(i, 'pivot');
      addTraversalStep({ description: `🎯 Finding minimum from position ${i} to ${n-1}`, highlightIndices: [i], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      
      updateElementState(minIdx, 'compared');
      addTraversalStep({ description: `📍 Current minimum candidate: ${arr[minIdx]} at index ${minIdx}`, highlightIndices: [minIdx, i], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      
      for (let j = i + 1; j < n; j++) {
        if (signal.aborted) return;
        while (isPaused) await sleep();
        
        updateElementState(j, 'active');
        localComparisons++;
        setComparisons(c => c + 1);
        addTraversalStep({ description: `🔍 Comparing ${arr[j]} with current minimum ${arr[minIdx]}`, highlightIndices: [j, minIdx], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
        
        await sleep();
        
        if (arr[j] < arr[minIdx]) {
          if (minIdx !== i) updateElementState(minIdx, 'default');
          minIdx = j;
          updateElementState(j, 'compared');
          addTraversalStep({ description: `🎉 New minimum found: ${arr[j]} at index ${j}`, highlightIndices: [j, i], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
        } else {
          addTraversalStep({ description: `❌ ${arr[j]} is not smaller than ${arr[minIdx]}`, highlightIndices: [j, minIdx], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
        }
        
        if (j !== minIdx && j !== i) updateElementState(j, 'default');
        await sleep();
      }
      
      if (minIdx !== i) {
        updateElementStates([i, minIdx], 'swap');
        addTraversalStep({ description: `🔄 Swapping ${arr[i]} with ${arr[minIdx]} (placing minimum in position ${i})`, highlightIndices: [i, minIdx], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
        
        [arr[i], arr[minIdx]] = [arr[minIdx], arr[i]];
        setArray([...arr]);
        localSwaps++;
        setSwaps(s => s + 1);
        
        await sleep();
      } else {
        addTraversalStep({ description: `✅ Element ${arr[i]} is already in correct position`, highlightIndices: [i], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      }
      
      updateElementState(i, 'ok');
      updateElementState(minIdx, 'default');
      addTraversalStep({ description: `🎯 Position ${i} is now sorted with value ${arr[i]}`, highlightIndices: [i], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      
      await sleep();
    }

    updateElementState(n - 1, 'ok');
    addTraversalStep({ description: "🎉 Selection Sort completed! All elements are now sorted.", arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
  };

  // Enhanced Insertion Sort implementation
  const insertionSort = async (signal) => {
    let arr = [...array];
    let n = arr.length;
    let localComparisons = 0, localSwaps = 0, localIterations = 0;
    addTraversalStep({ description: "🚀 Starting Insertion Sort - building sorted array one element at a time", arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });

    updateElementState(0, 'ok');
    addTraversalStep({ description: `✅ First element ${arr[0]} is considered sorted`, highlightIndices: [0], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });

    for (let i = 1; i < n; i++) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(s => s + 1);
      let key = arr[i];
      let j = i - 1;
      
      updateElementState(i, 'active');
      addTraversalStep({ description: `📌 Inserting element ${key} from position ${i} into sorted portion`, highlightIndices: [i], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      
      await sleep();
      
      while (j >= 0 && arr[j] > key) {
        if (signal.aborted) return;
        while (isPaused) await sleep();
        
        localComparisons++;
        setComparisons(c => c + 1);
        updateElementStates([j, j + 1], 'swap');
        addTraversalStep({ description: `🔍 ${arr[j]} > ${key}, shifting ${arr[j]} to the right`, highlightIndices: [j, j + 1], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
        
        arr[j + 1] = arr[j];
        setArray([...arr]);
        localSwaps++;
        setSwaps(s => s + 1);
        
        await sleep();
        j = j - 1;
      }
      
      arr[j + 1] = key;
      setArray([...arr]);
      
      updateElementState(j + 1, 'found');
      addTraversalStep({ description: `🎯 Inserted ${key} at position ${j + 1}`, highlightIndices: [j + 1], arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      
      await sleep();
      
      for (let k = 0; k <= i; k++) {
        updateElementState(k, 'ok');
      }
      addTraversalStep({ description: `✅ Elements 0 to ${i} are now sorted`, highlightIndices: Array.from({length: i + 1}, (_, k) => k), arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
      
      await sleep();
    }

    addTraversalStep({ description: "🎉 Insertion Sort completed! All elements are now sorted.", arrayState: arr, comparisons: localComparisons, swaps: localSwaps, iterations: localIterations });
  };

  // Enhanced Merge Sort implementation
  const mergeSort = async (signal) => {
    let arr = [...array];
    let stats = { comparisons: 0, swaps: 0, iterations: 0 };
    addTraversalStep({ description: "🚀 Starting Merge Sort - divide and conquer approach", arrayState: arr, ...stats });
    addTraversalStep({ description: `📊 Initial array: [${arr.join(', ')}]`, arrayState: arr, ...stats });

    await mergeSortHelper(arr, 0, arr.length - 1, signal, stats);

    addTraversalStep({ description: "🎉 Merge Sort completed! All elements are now sorted.", arrayState: arr, ...stats });
  };

  // Recursive helper function for merge sort
  const mergeSortHelper = async (arr, left, right, signal, stats) => {
    if (left >= right) return;
    const mid = Math.floor((left + right) / 2);

    addTraversalStep({ description: `✂️ Splitting array from index ${left} to ${right} at midpoint ${mid}`, arrayState: arr, ...stats });

    for (let i = left; i <= right; i++) {
      if (i <= mid) {
        updateElementState(i, 'active');
      } else {
        updateElementState(i, 'compared');
      }
    }
    await sleep();

    await mergeSortHelper(arr, left, mid, signal, stats);
    await mergeSortHelper(arr, mid + 1, right, signal, stats);

    await merge(arr, left, mid, right, signal, stats);
  };

  // Merge function for merge sort
  const merge = async (arr, left, mid, right, signal, stats) => {
    addTraversalStep({ description: `🔄 Merging sorted subarrays from ${left}-${mid} and ${mid+1}-${right}`, arrayState: arr, ...stats });

    let temp = [];
    let i = left;
    let j = mid + 1;

    while (i <= mid && j <= right) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      stats.comparisons++;
      setComparisons(c => c + 1);
      updateElementStates([i, j], 'active');
      addTraversalStep({ description: `🔍 Comparing ${arr[i]} (index ${i}) and ${arr[j]} (index ${j})`, highlightIndices: [i, j], arrayState: arr, ...stats });
      
      await sleep();
      
      if (arr[i] <= arr[j]) {
        temp.push(arr[i]);
        i++;
      } else {
        temp.push(arr[j]);
        j++;
      }
    }

    while (i <= mid) {
      temp.push(arr[i]);
      i++;
    }

    while (j <= right) {
      temp.push(arr[j]);
      j++;
    }

    for (let k = left; k <= right; k++) {
      arr[k] = temp[k - left];
      setArray([...arr]);
      updateElementState(k, 'swap');
      addTraversalStep({ description: `📝 Writing ${arr[k]} to position ${k}`, highlightIndices: [k], arrayState: arr, ...stats });
      await sleep();
    }

    for (let k = left; k <= right; k++) {
      updateElementState(k, 'ok');
    }
    addTraversalStep({ description: `✅ Segment ${left}-${right} is now sorted`, highlightIndices: Array.from({length: right - left + 1}, (_, k) => k + left), arrayState: arr, ...stats });
  };

  // Enhanced Quick Sort implementation
  const quickSort = async (signal) => {
    let arr = [...array];
    let stats = { comparisons: 0, swaps: 0, iterations: 0 };
    addTraversalStep({ description: "🚀 Starting Quick Sort - divide and conquer with partitioning", arrayState: arr, ...stats });
    addTraversalStep({ description: `📊 Initial array: [${arr.join(', ')}]`, arrayState: arr, ...stats });

    await quickSortHelper(arr, 0, arr.length - 1, signal, stats);

    for(let i = 0; i < arr.length; i++) updateElementState(i, 'ok');
    addTraversalStep({ description: "🎉 Quick Sort completed! All elements are now sorted.", arrayState: arr, ...stats });
  };

  // Recursive helper function for quick sort
  const quickSortHelper = async (arr, low, high, signal, stats) => {
    if (low < high) {
      if (signal.aborted) return;
      let pi = await partition(arr, low, high, signal, stats);
      
      await quickSortHelper(arr, low, pi - 1, signal, stats);
      await quickSortHelper(arr, pi + 1, high, signal, stats);
    } else if (low >= 0 && low < arr.length) {
      updateElementState(low, 'ok');
    }
  };

  // Partition function for quick sort
  const partition = async (arr, low, high, signal, stats) => {
    stats.iterations++;
    setIterations(s => s + 1);
    let pivot = arr[high];

    updateElementState(high, 'pivot');
    addTraversalStep({ description: `🎯 Selecting pivot: ${pivot} at index ${high}`, highlightIndices: [high], arrayState: arr, ...stats });
    await sleep();

    let i = low - 1;

    for (let j = low; j < high; j++) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      stats.comparisons++;
      setComparisons(c => c + 1);
      updateElementStates([j, high], 'active');
      addTraversalStep({ description: `🔍 Comparing ${arr[j]} with pivot ${pivot}`, highlightIndices: [j, high], arrayState: arr, ...stats });
      
      await sleep();
      
      if (arr[j] < pivot) {
        i++;
        
        if (i !== j) {
          updateElementStates([i, j], 'swap');
          addTraversalStep({ description: `🔄 Swapping ${arr[i]} and ${arr[j]} (${arr[j]} < ${pivot})`, highlightIndices: [i, j], arrayState: arr, ...stats });
          
          [arr[i], arr[j]] = [arr[j], arr[i]];
          setArray([...arr]);
          stats.swaps++;
          setSwaps(s => s + 1);
          
          await sleep();
        }
      } else {
        addTraversalStep({ description: `❌ ${arr[j]} >= ${pivot}, no swap needed`, highlightIndices: [j, high], arrayState: arr, ...stats });
      }
      
      updateElementStates([j, high], 'default');
      updateElementState(high, 'pivot'); // keep pivot highlighted
    }

    i++;
    if (i !== high) {
      updateElementStates([i, high], 'swap');
      addTraversalStep({ description: `🔄 Placing pivot in correct position: swapping ${arr[i]} and ${arr[high]}`, highlightIndices: [i, high], arrayState: arr, ...stats });
      
      [arr[i], arr[high]] = [arr[high], arr[i]];
      setArray([...arr]);
      stats.swaps++;
      setSwaps(s => s + 1);
      
      await sleep();
    } else {
      addTraversalStep({ description: `✅ Pivot ${pivot} is already in correct position at index ${i}`, highlightIndices: [i], arrayState: arr, ...stats });
    }

    updateElementState(i, 'ok');
    addTraversalStep({ description: `🎯 Pivot ${arr[i]} is now in its final position at index ${i}`, highlightIndices: [i], arrayState: arr, ...stats });

    return i;
  };

  // Enhanced Heap Sort implementation
  const heapSort = async (signal) => {
    let arr = [...array];
    let n = arr.length;
    let stats = { comparisons: 0, swaps: 0, iterations: 0 };
    addTraversalStep({ description: "🚀 Starting Heap Sort - building a max heap", arrayState: arr, ...stats });
    addTraversalStep({ description: `📊 Initial array: [${arr.join(', ')}]`, arrayState: arr, ...stats });

    for (let i = Math.floor(n / 2) - 1; i >= 0; i--) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      addTraversalStep({ description: `🏗️ Building max heap: heapifying node at index ${i}`, arrayState: arr, ...stats });
      await heapify(arr, n, i, signal, stats);
    }

    for (let i = n - 1; i > 0; i--) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      stats.iterations++;
      setIterations(s => s + 1);
      
      updateElementStates([0, i], 'swap');
      addTraversalStep({ description: `🔄 Moving root ${arr[0]} to end at position ${i}`, highlightIndices: [0, i], arrayState: arr, ...stats });
      
      [arr[0], arr[i]] = [arr[i], arr[0]];
      setArray([...arr]);
      stats.swaps++;
      setSwaps(s => s + 1);
      
      await sleep();
      
      updateElementState(i, 'ok');
      addTraversalStep({ description: `✅ Element ${arr[i]} is now in its final position`, highlightIndices: [i], arrayState: arr, ...stats });
      
      addTraversalStep({ description: `⚙️ Heapifying the reduced heap of size ${i}`, arrayState: arr, ...stats });
      await heapify(arr, i, 0, signal, stats);
    }

    updateElementState(0, 'ok');
    addTraversalStep({ description: "🎉 Heap Sort completed! All elements are now sorted.", arrayState: arr, ...stats });
  };

  // Heapify function for heap sort
  const heapify = async (arr, n, i, signal, stats) => {
    let largest = i;
    let left = 2 * i + 1;
    let right = 2 * i + 2;
    let indicesToHighlight = [i];
    if (left < n) indicesToHighlight.push(left);
    if (right < n) indicesToHighlight.push(right);

    updateElementStates(indicesToHighlight, 'active');
    addTraversalStep({ description: `🔍 Comparing node ${arr[i]} at index ${i} with its children`, highlightIndices: indicesToHighlight, arrayState: arr, ...stats });

    await sleep();

    if (left < n) {
      stats.comparisons++;
      setComparisons(c => c+1);
      if (arr[left] > arr[largest]) {
        largest = left;
        addTraversalStep({ description: `📈 Left child ${arr[left]} is larger than current largest ${arr[i]}`, highlightIndices: [left, i], arrayState: arr, ...stats });
      }
    }

    if (right < n) {
      stats.comparisons++;
      setComparisons(c => c+1);
      if (arr[right] > arr[largest]) {
        largest = right;
        addTraversalStep({ description: `📈 Right child ${arr[right]} is larger than current largest ${arr[i]}`, highlightIndices: [right, i], arrayState: arr, ...stats });
      }
    }

    if (largest !== i) {
      updateElementStates([i, largest], 'swap');
      addTraversalStep({ description: `🔄 Swapping ${arr[i]} and ${arr[largest]}`, highlightIndices: [i, largest], arrayState: arr, ...stats });
      
      [arr[i], arr[largest]] = [arr[largest], arr[i]];
      setArray([...arr]);
      stats.swaps++;
      setSwaps(s => s + 1);
      
      await sleep();
      
      addTraversalStep({ description: `⚙️ Recursively heapifying sub-tree rooted at index ${largest}`, arrayState: arr, ...stats });
      await heapify(arr, n, largest, signal, stats);
    } else {
      addTraversalStep({ description: `✅ Node ${arr[i]} is already the largest in its sub-tree`, highlightIndices: [i], arrayState: arr, ...stats });
    }

    updateElementStates(indicesToHighlight, 'default');
  };

  // Enhanced Linear Search implementation
  const linearSearch = async (signal) => {
    let arr = [...array];
    const targetVal = parseInt(target) || arr[Math.floor(Math.random() * arr.length)];
    setTarget(targetVal.toString());
    setSearchResult(null);
    resetElementStates();
    let localComparisons = 0, localIterations = 0;

    addTraversalStep({ description: `🔍 Starting Linear Search for target value: ${targetVal}`, arrayState: arr, comparisons: localComparisons, swaps: 0, iterations: localIterations });

    for (let i = 0; i < arr.length; i++) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(s => s + 1);
      updateElementState(i, 'active');
      localComparisons++;
      setComparisons(c => c + 1);
      addTraversalStep({ description: `🔍 Checking position ${i}: comparing ${arr[i]} with target ${targetVal}`, highlightIndices: [i], arrayState: arr, comparisons: localComparisons, swaps: 0, iterations: localIterations });
      
      await sleep();
      
      if (arr[i] === targetVal) {
        updateElementState(i, 'found');
        addTraversalStep({ description: `🎉 Target ${targetVal} found at index ${i}!`, highlightIndices: [i], arrayState: arr, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        setSearchResult(i);
        return;
      }
      
      updateElementState(i, 'visited');
      addTraversalStep({ description: `❌ ${arr[i]} ≠ ${targetVal}, continuing search...`, highlightIndices: [i], arrayState: arr, comparisons: localComparisons, swaps: 0, iterations: localIterations });
    }

    addTraversalStep({ description: `😞 Target ${targetVal} not found in the array`, arrayState: arr, comparisons: localComparisons, swaps: 0, iterations: localIterations });
    setSearchResult(-1);
  };

  // Enhanced Binary Search implementation
  const binarySearch = async (signal) => {
    const sortedArray = [...array].sort((a, b) => a - b);
    setArray(sortedArray);
    resetElementStates();
    const targetVal = parseInt(target) || sortedArray[Math.floor(Math.random() * sortedArray.length)];
    setTarget(targetVal.toString());
    setSearchResult(null);

    let left = 0;
    let right = sortedArray.length - 1;
    let localComparisons = 0, localIterations = 0;

    addTraversalStep({ description: `🚀 Starting Binary Search for target: ${targetVal} (array is sorted)`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });

    while (left <= right) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(s => s + 1);
      const mid = Math.floor((left + right) / 2);
      
      for (let i = 0; i < sortedArray.length; i++) {
        if (i < left || i > right) updateElementState(i, 'visited');
        else if (i === mid) updateElementState(i, 'active');
        else updateElementState(i, 'compared');
      }
      
      localComparisons++;
      setComparisons(c => c + 1);
      addTraversalStep({ description: `🎯 Search range [${left}, ${right}], checking middle index ${mid} (value: ${sortedArray[mid]})`, highlightIndices: [mid], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
      
      await sleep();
      
      if (sortedArray[mid] === targetVal) {
        updateElementState(mid, 'found');
        addTraversalStep({ description: `🎉 Target ${targetVal} found at index ${mid}!`, highlightIndices: [mid], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        setSearchResult(mid);
        return;
      } else if (sortedArray[mid] < targetVal) {
        addTraversalStep({ description: `📈 ${sortedArray[mid]} < ${targetVal}, searching right half [${mid + 1}, ${right}]`, highlightIndices: [mid], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        left = mid + 1;
      } else {
        addTraversalStep({ description: `📉 ${sortedArray[mid]} > ${targetVal}, searching left half [${left}, ${mid - 1}]`, highlightIndices: [mid], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        right = mid - 1;
      }
      
      await sleep();
    }

    addTraversalStep({ description: `😞 Target ${targetVal} not found in the array`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
    setSearchResult(-1);
  };

  // Enhanced Jump Search implementation
  const jumpSearch = async (signal) => {
    const sortedArray = [...array].sort((a, b) => a - b);
    setArray(sortedArray);
    resetElementStates();
    const targetVal = parseInt(target) || sortedArray[Math.floor(Math.random() * sortedArray.length)];
    setTarget(targetVal.toString());
    setSearchResult(null);

    const n = sortedArray.length;
    let step = Math.floor(Math.sqrt(n));
    let prev = 0;
    let localComparisons = 0, localIterations = 0;

    addTraversalStep({ description: `🚀 Starting Jump Search for target: ${targetVal} (array is sorted)`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
    addTraversalStep({ description: `📏 Jump step size: ${step} (sqrt(${n}) ≈ ${step})`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });

    for (let i = 0; i < Math.min(step, n); i++) updateElementState(i, 'compared');
    addTraversalStep({ description: `🔍 Checking block [0, ${Math.min(step, n)-1}]`, highlightIndices: Array.from({length: Math.min(step, n)}, (_, i) => i), arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });

    while (sortedArray[Math.min(step, n) - 1] < targetVal) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(i => i + 1);
      localComparisons++;
      setComparisons(c => c + 1);
      prev = step;
      step += Math.floor(Math.sqrt(n));
      
      if (prev >= n) {
        addTraversalStep({ description: `❌ Previous block (${prev}) exceeds array length, target not found`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        setSearchResult(-1);
        return;
      }
      
      for (let i = prev; i < Math.min(step, n); i++) updateElementState(i, 'compared');
      addTraversalStep({ description: `🔍 Checking block [${prev}, ${Math.min(step, n)-1}]`, highlightIndices: Array.from({length: Math.min(step, n) - prev}, (_, i) => i + prev), arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
      
      await sleep();
    }

    addTraversalStep({ description: `🔎 Performing linear search within block [${prev}, ${Math.min(step, n)-1}]`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });

    while (sortedArray[prev] < targetVal) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(i => i + 1);
      localComparisons++;
      setComparisons(c => c + 1);
      
      updateElementState(prev, 'active');
      addTraversalStep({ description: `🔍 Checking element at index ${prev}: ${sortedArray[prev]} < ${targetVal}`, highlightIndices: [prev], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
      
      await sleep();
      
      prev++;
      
      if (prev === Math.min(step, n)) {
        addTraversalStep({ description: `❌ Reached end of block, target not found`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        setSearchResult(-1);
        return;
      }
    }

    localComparisons++;
    setComparisons(c => c + 1);
    updateElementState(prev, 'active');
    addTraversalStep({ description: `🔍 Checking element at index ${prev}: ${sortedArray[prev]} == ${targetVal}?`, highlightIndices: [prev], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });

    await sleep();

    if (sortedArray[prev] === targetVal) {
      updateElementState(prev, 'found');
      addTraversalStep({ description: `🎉 Target ${targetVal} found at index ${prev}!`, highlightIndices: [prev], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
      setSearchResult(prev);
      return;
    }

    addTraversalStep({ description: `❌ Element at index ${prev} is ${sortedArray[prev]} ≠ ${targetVal}`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
    setSearchResult(-1);
  };

  // Enhanced Interpolation Search implementation
  const interpolationSearch = async (signal) => {
    const sortedArray = [...array].sort((a, b) => a - b);
    setArray(sortedArray);
    resetElementStates();
    const targetVal = parseInt(target) || sortedArray[Math.floor(Math.random() * sortedArray.length)];
    setTarget(targetVal.toString());
    setSearchResult(null);

    let low = 0;
    let high = sortedArray.length - 1;
    let localComparisons = 0, localIterations = 0;

    addTraversalStep({ description: `🚀 Starting Interpolation Search for target: ${targetVal} (array is sorted)`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });

    while (low <= high && targetVal >= sortedArray[low] && targetVal <= sortedArray[high]) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(i => i + 1);
      
      if (low === high) {
        localComparisons++;
        setComparisons(c => c + 1);
        updateElementState(low, 'active');
        addTraversalStep({ description: `🔍 Only one element left: checking index ${low} (${sortedArray[low]} == ${targetVal}?)`, highlightIndices: [low], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        
        await sleep();
        
        if (sortedArray[low] === targetVal) {
          updateElementState(low, 'found');
          addTraversalStep({ description: `🎉 Target ${targetVal} found at index ${low}!`, highlightIndices: [low], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
          setSearchResult(low);
          return;
        }
        
        addTraversalStep({ description: `❌ Single element ${sortedArray[low]} ≠ ${targetVal}`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        setSearchResult(-1);
        return;
      }
      
      const pos = low + Math.floor(((targetVal - sortedArray[low]) * (high - low)) / (sortedArray[high] - sortedArray[low]));
      
      for (let i = low; i <= high; i++) {
        if (i === pos) updateElementState(i, 'active');
        else updateElementState(i, 'compared');
      }
      
      localComparisons++;
      setComparisons(c => c + 1);
      addTraversalStep({ description: `🧮 Calculated position: ${pos} (using interpolation formula)`, highlightIndices: [pos], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
      addTraversalStep({ description: `🔍 Checking element at index ${pos}: ${sortedArray[pos]} == ${targetVal}?`, highlightIndices: [pos], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
      
      await sleep();
      
      if (sortedArray[pos] === targetVal) {
        updateElementState(pos, 'found');
        addTraversalStep({ description: `🎉 Target ${targetVal} found at index ${pos}!`, highlightIndices: [pos], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        setSearchResult(pos);
        return;
      }
      
      if (sortedArray[pos] < targetVal) {
        addTraversalStep({ description: `📈 ${sortedArray[pos]} < ${targetVal}, searching right half [${pos + 1}, ${high}]`, highlightIndices: [pos], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        low = pos + 1;
      } else {
        addTraversalStep({ description: `📉 ${sortedArray[pos]} > ${targetVal}, searching left half [${low}, ${pos - 1}]`, highlightIndices: [pos], arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
        high = pos - 1;
      }
      
      await sleep();
    }

    addTraversalStep({ description: `😞 Target ${targetVal} not found in the array`, arrayState: sortedArray, comparisons: localComparisons, swaps: 0, iterations: localIterations });
    setSearchResult(-1);
  };

  // Graph BFS implementation
  const graphBFS = async (signal) => {
    if (!startNode || !graph[startNode]) {
      addTraversalStep({ description: "❌ Please select a valid start node" });
      return;
    }
    let localComparisons = 0, localIterations = 0;
    addTraversalStep({ description: `🚀 Starting Graph BFS from node ${startNode}`, highlightIndices: [startNode], comparisons: localComparisons, iterations: localIterations });

    let visited = new Set([startNode]);
    let queue = [startNode];

    while (queue.length > 0) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(i => i + 1);
      let currentNode = queue.shift();
      
      addTraversalStep({ description: `🔍 Processing node ${currentNode}`, highlightIndices: [currentNode], comparisons: localComparisons, iterations: localIterations });
      await sleep();
      
      for (let neighbor in graph[currentNode]) {
        if (signal.aborted) return;
        while (isPaused) await sleep();
        
        localComparisons++;
        setComparisons(c => c + 1);
        
        if (!visited.has(neighbor)) {
          visited.add(neighbor);
          queue.push(neighbor);
          addTraversalStep({ description: `📌 Discovered node ${neighbor} from ${currentNode}`, highlightIndices: [currentNode, neighbor], comparisons: localComparisons, iterations: localIterations });
          await sleep();
        } else {
          addTraversalStep({ description: `⏩ Already visited node ${neighbor} from ${currentNode}`, highlightIndices: [currentNode, neighbor], comparisons: localComparisons, iterations: localIterations });
        }
      }
    }

    addTraversalStep({ description: "🎉 Graph BFS completed!", comparisons: localComparisons, iterations: localIterations });
  };

  // Graph DFS implementation
  const graphDFS = async (signal) => {
    if (!startNode || !graph[startNode]) {
      addTraversalStep({ description: "❌ Please select a valid start node" });
      return;
    }
    let localComparisons = 0, localIterations = 0;
    addTraversalStep({ description: `🚀 Starting Graph DFS from node ${startNode}`, highlightIndices: [startNode], comparisons: localComparisons, iterations: localIterations });

    let visited = new Set();
    let stack = [startNode];
    
    while (stack.length > 0) {
      if (signal.aborted) return;
      while (isPaused) await sleep();

      let currentNode = stack.pop();

      if (!visited.has(currentNode)) {
        visited.add(currentNode);
        localIterations++;
        setIterations(i => i + 1);

        addTraversalStep({ description: `🔍 Processing node ${currentNode}`, highlightIndices: [currentNode], comparisons: localComparisons, iterations: localIterations });
        await sleep();

        const neighbors = Object.keys(graph[currentNode] || {}).reverse();
        for (let neighbor of neighbors) {
          if (signal.aborted) return;
          localComparisons++;
          setComparisons(c => c + 1);
          if (!visited.has(neighbor)) {
            stack.push(neighbor);
            addTraversalStep({ description: `📌 Adding undiscovered node ${neighbor} to the stack`, highlightIndices: [currentNode, neighbor], comparisons: localComparisons, iterations: localIterations });
          } else {
            addTraversalStep({ description: `⏩ Already visited node ${neighbor} from ${currentNode}`, highlightIndices: [currentNode, neighbor], comparisons: localComparisons, iterations: localIterations });
          }
        }
      }
    }
    addTraversalStep({ description: "🎉 Graph DFS completed!", comparisons: localComparisons, iterations: localIterations });
  };

  // Dijkstra's algorithm implementation
  const graphDijkstra = async (signal) => {
    if (!startNode || !graph[startNode]) {
      addTraversalStep({ description: "❌ Please select a valid start node"});
      return;
    }
    let localComparisons = 0, localIterations = 0;
    addTraversalStep({ description: `🚀 Starting Dijkstra's Algorithm from node ${startNode}`, highlightIndices: [startNode], comparisons: localComparisons, iterations: localIterations });

    let distances = {};
    let visited = new Set();

    for (let node in graph) {
      distances[node] = Infinity;
    }
    distances[startNode] = 0;

    addTraversalStep({ description: `📊 Initialized distances`, graphState: { distances }, comparisons: localComparisons, iterations: localIterations });

    while (true) {
      if (signal.aborted) return;
      while (isPaused) await sleep();
      
      localIterations++;
      setIterations(i => i + 1);
      
      let minNode = null;
      for (let node in graph) {
        if (!visited.has(node) && (minNode === null || distances[node] < distances[minNode])) {
          minNode = node;
        }
      }
      
      if (minNode === null) break;
      
      visited.add(minNode);
      addTraversalStep({ description: `🔍 Visiting node ${minNode} with distance ${distances[minNode]}`, highlightIndices: [minNode], graphState: { distances }, comparisons: localComparisons, iterations: localIterations });
      await sleep();
      
      for (let neighbor in graph[minNode]) {
        if (signal.aborted) return;
        while (isPaused) await sleep();
        
        localComparisons++;
        setComparisons(c => c + 1);
        
        let newDistance = distances[minNode] + graph[minNode][neighbor];
        if (newDistance < distances[neighbor]) {
          distances[neighbor] = newDistance;
          addTraversalStep({ description: `📈 Updated distance to ${neighbor}: ${newDistance} (via ${minNode})`, highlightIndices: [minNode, neighbor], graphState: { distances }, comparisons: localComparisons, iterations: localIterations });
          await sleep();
        } else {
          addTraversalStep({ description: `📉 No update for ${neighbor}: ${distances[neighbor]} ≤ ${newDistance}`, highlightIndices: [minNode, neighbor], graphState: { distances }, comparisons: localComparisons, iterations: localIterations });
        }
      }
    }

    addTraversalStep({ description: `🎉 Dijkstra's Algorithm completed!`, graphState: { distances, final: true }, comparisons: localComparisons, iterations: localIterations });
  };

  // Bellman-Ford algorithm implementation
  const bellmanFord = async (signal) => {
    if (!startNode || !graph[startNode]) {
      addTraversalStep({ description: "❌ Please select a valid start node" });
      return;
    }
    let localComparisons = 0, localIterations = 0;
    addTraversalStep({ description: `🚀 Starting Bellman-Ford Algorithm from node ${startNode}`, highlightIndices: [startNode], comparisons: localComparisons, iterations: localIterations });

    let distances = {};
    for (let vertex in graph) {
      distances[vertex] = Infinity;
    }
    distances[startNode] = 0;

    addTraversalStep({ description: `📊 Initialized distances`, graphState: { distances }, comparisons: localComparisons, iterations: localIterations });

    const V = Object.keys(graph).length;
    for (let i = 0; i < V - 1; i++) {
      if (signal.aborted) return;
      while (isPaused) await sleep();

      localIterations++;
      setIterations(it => it + 1);
      addTraversalStep({ description: `🔄 Iteration ${i + 1}/${V - 1}`, graphState: { distances }, comparisons: localComparisons, iterations: localIterations });

      for (let u in graph) {
        for (let v in graph[u]) {
          if (signal.aborted) return;
          while (isPaused) await sleep();

          localComparisons++;
          setComparisons(c => c + 1);

          if (distances[u] !== Infinity && distances[u] + graph[u][v] < distances[v]) {
            distances[v] = distances[u] + graph[u][v];
            addTraversalStep({ description: `📈 Updated distance to ${v}: ${distances[v]} (via ${u})`, highlightIndices: [u, v], graphState: { distances }, comparisons: localComparisons, iterations: localIterations });
            await sleep();
          }
        }
      }
    }

    // Check for negative weight cycles
    for (let u in graph) {
      for (let v in graph[u]) {
        if (distances[u] !== Infinity && distances[u] + graph[u][v] < distances[v]) {
          addTraversalStep({ description: `🚨 Negative weight cycle detected!`, highlightIndices: [u, v], graphState: { distances }, comparisons: localComparisons, iterations: localIterations });
          return;
        }
      }
    }

    addTraversalStep({ description: `🎉 Bellman-Ford Algorithm completed!`, graphState: { distances, final: true }, comparisons: localComparisons, iterations: localIterations });
  };

  // Prim's algorithm implementation
  const prims = async (signal) => {
    if (Object.keys(graph).length === 0) {
      addTraversalStep({ description: "❌ Graph is empty" });
      return;
    }
    let localComparisons = 0, localIterations = 0;
    addTraversalStep({ description: `🚀 Starting Prim's Algorithm`, comparisons: localComparisons, iterations: localIterations });

    let mst = [];
    let visited = new Set();
    let startNode = Object.keys(graph)[0];
    visited.add(startNode);

    while (visited.size < Object.keys(graph).length) {
      if (signal.aborted) return;
      while (isPaused) await sleep();

      localIterations++;
      setIterations(it => it + 1);

      let minEdge = null;
      for (let u of visited) {
        for (let v in graph[u]) {
          if (!visited.has(v)) {
            localComparisons++;
            setComparisons(c => c + 1);
            if (minEdge === null || graph[u][v] < minEdge.weight) {
              minEdge = { u, v, weight: graph[u][v] };
            }
          }
        }
      }

      if (minEdge) {
        mst.push(minEdge);
        visited.add(minEdge.v);
        addTraversalStep({ description: `➕ Adding edge ${minEdge.u}-${minEdge.v} to MST`, highlightIndices: [minEdge.u, minEdge.v], graphState: { mst }, comparisons: localComparisons, iterations: localIterations });
        await sleep();
      } else {
        break; // No more edges to add
      }
    }

    addTraversalStep({ description: `🎉 Prim's Algorithm completed!`, graphState: { mst, final: true }, comparisons: localComparisons, iterations: localIterations });
  };

  // Kruskal's algorithm implementation
  const kruskal = async (signal) => {
    if (Object.keys(graph).length === 0) {
      addTraversalStep({ description: "❌ Graph is empty" });
      return;
    }
    let localComparisons = 0, localIterations = 0;
    addTraversalStep({ description: `🚀 Starting Kruskal's Algorithm`, comparisons: localComparisons, iterations: localIterations });

    let edges = [];
    for (let u in graph) {
      for (let v in graph[u]) {
        edges.push({ u, v, weight: graph[u][v] });
      }
    }
    edges.sort((a, b) => a.weight - b.weight);

    let parent = {};
    Object.keys(graph).forEach(node => parent[node] = node);

    function find(i) {
      if (parent[i] === i) return i;
      return find(parent[i]);
    }

    function union(i, j) {
      let rootI = find(i);
      let rootJ = find(j);
      if (rootI !== rootJ) {
        parent[rootJ] = rootI;
        return true;
      }
      return false;
    }

    let mst = [];
    for (let edge of edges) {
      if (signal.aborted) return;
      while (isPaused) await sleep();

      localIterations++;
      setIterations(it => it + 1);
      localComparisons++;
      setComparisons(c => c + 1);

      if (union(edge.u, edge.v)) {
        mst.push(edge);
        addTraversalStep({ description: `➕ Adding edge ${edge.u}-${edge.v} to MST`, highlightIndices: [edge.u, edge.v], graphState: { mst }, comparisons: localComparisons, iterations: localIterations });
        await sleep();
      } else {
        addTraversalStep({ description: `➖ Discarding edge ${edge.u}-${edge.v} (creates a cycle)`, highlightIndices: [edge.u, edge.v], graphState: { mst }, comparisons: localComparisons, iterations: localIterations });
        await sleep();
      }
    }

    addTraversalStep({ description: `🎉 Kruskal's Algorithm completed!`, graphState: { mst, final: true }, comparisons: localComparisons, iterations: localIterations });
  };

  // Run the selected sorting algorithm
  const runSortingAlgorithm = async (name, signal) => {
    switch (name) {
      case 'bubble':
        await bubbleSort(signal);
        break;
      case 'selection':
        await selectionSort(signal);
        break;
      case 'insertion':
        await insertionSort(signal);
        break;
      case 'merge':
        await mergeSort(signal);
        break;
      case 'quick':
        await quickSort(signal);
        break;
      case 'heap':
        await heapSort(signal);
        break;
      default:
        addTraversalStep({ description: `⚠️ Algorithm ${name} not yet implemented`, arrayState: array, comparisons, swaps, iterations });
    }
  };

  // Run the selected search algorithm
  const runSearchAlgorithm = async (name, signal) => {
    switch (name) {
      case 'linear':
        await linearSearch(signal);
        break;
      case 'binary':
        await binarySearch(signal);
        break;
      case 'jump':
        await jumpSearch(signal);
        break;
      case 'interpolation':
        await interpolationSearch(signal);
        break;
      default:
        addTraversalStep({ description: `⚠️ Algorithm ${name} not yet implemented`, arrayState: array, comparisons, swaps, iterations });
    }
  };

  // Run the selected graph algorithm
  const runGraphAlgorithm = async (name, signal) => {
    switch (name) {
      case 'bfs':
        await graphBFS(signal);
        break;
      case 'dfs':
        await graphDFS(signal);
        break;
      case 'dijkstra':
        await graphDijkstra(signal);
        break;
      case 'bellmanford':
        await bellmanFord(signal);
        break;
      case 'prims':
        await prims(signal);
        break;
      case 'kruskal':
        await kruskal(signal);
        break;
      default:
        addTraversalStep({ description: `⚠️ Algorithm ${name} not yet implemented` });
    }
  };

  // Render array elements with enhanced styling
  const renderArrayElements = () => {
    return arrayElements.map((element, index) => {
      let bgClass = 'bg-gradient-to-br from-green-400 to-green-600';
      let transform = 'scale-100';
      let shadow = 'shadow-md';
      if (element.state === 'active') {
        bgClass = 'bg-gradient-to-br from-blue-400 to-blue-600';
        transform = 'scale-110';
        shadow = 'shadow-lg';
      } else if (element.state === 'compared') {
        bgClass = 'bg-gradient-to-br from-yellow-400 to-yellow-600';
        transform = 'scale-105';
        shadow = 'shadow-lg';
      } else if (element.state === 'swap') {
        bgClass = 'bg-gradient-to-br from-pink-400 to-pink-600';
        transform = 'scale-110';
        shadow = 'shadow-xl';
      } else if (element.state === 'ok') {
        bgClass = 'bg-gradient-to-br from-teal-400 to-teal-600';
        shadow = 'shadow-lg';
      } else if (element.state === 'found') {
        bgClass = 'bg-gradient-to-br from-purple-400 to-purple-600';
        transform = 'scale-115';
        shadow = 'shadow-2xl';
      } else if (element.state === 'visited') {
        bgClass = 'bg-gradient-to-br from-gray-400 to-gray-600';
        shadow = 'shadow-sm';
      } else if (element.state === 'pivot') {
        bgClass = 'bg-gradient-to-br from-orange-400 to-orange-600';
        transform = 'scale-105';
        shadow = 'shadow-lg';
      }
      
      return (
        <div
          key={index}
          className={`w-14 flex items-center justify-center rounded-lg text-white font-bold text-lg transition-all duration-300 ${bgClass} ${shadow} transform ${transform}`}
          style={{ 
            height: `${Math.max(element.value * 3, 40)}px`,
            minHeight: '40px'
          }}
          title={`Value: ${element.value}, Index: ${index}, State: ${element.state}`}
        >
          {element.value}
        </div>
      );
    });
  };

  // Render the active visualization based on algorithm type
  const renderActiveVisualization = () => {
    const [group] = selectedAlgorithm.split(':');
    if (group === 'sort' || group === 'search') {
      return (
        <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl p-6 min-h-[500px] border border-gray-700 shadow-inner flex flex-col justify-between">
          <div>
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-bold text-blue-400">📊 Array Visualization</h3>
              <div className="text-sm text-gray-400">
                Elements: {array.length} | Max: {Math.max(...array, 0)} | Min: {Math.min(...array, 0)}
              </div>
            </div>
            <div className="flex flex-wrap gap-3 justify-center items-end mb-4" ref={arrayContainerRef}>
              {renderArrayElements()}
            </div>
          </div>
          <div>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-2 mt-4 text-xs">
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-blue-400 to-blue-600 rounded"></div><span>Active</span></div>
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-yellow-400 to-yellow-600 rounded"></div><span>Comparing</span></div>
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-pink-400 to-pink-600 rounded"></div><span>Swapping</span></div>
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-teal-400 to-teal-600 rounded"></div><span>Sorted</span></div>
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-purple-400 to-purple-600 rounded"></div><span>Found</span></div>
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-gray-400 to-gray-600 rounded"></div><span>Visited</span></div>
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-orange-400 to-orange-600 rounded"></div><span>Pivot</span></div>
              <div className="flex items-center space-x-2"><div className="w-3 h-3 bg-gradient-to-br from-green-400 to-green-600 rounded"></div><span>Default</span></div>
            </div>
            {searchResult !== null && (
              <div className="mt-4 p-4 bg-gradient-to-r from-gray-800 to-gray-750 rounded-lg border border-gray-600 text-center">
                {searchResult >= 0 ? <p className="text-purple-400 font-bold text-lg">🎉 Target found at index {searchResult}!</p> : <p className="text-red-400 font-bold text-lg">😞 Target not found in the array</p>}
              </div>
            )}
          </div>
        </div>
      );
    }

    if (group === 'graph') {
      return <GraphCanvas 
                nodes={graphNodes}
                edges={graphEdges}
                positions={nodePositions}
                setPositions={setNodePositions}
                currentStepData={traversalSteps[currentStep]}
             />;
    }

    return null;
  };

  // GraphCanvas component for visualizing graph algorithms
  const GraphCanvas = ({ nodes, edges, positions, setPositions, currentStepData }) => {
    const containerRef = useRef(null);
    const [draggingNode, setDraggingNode] = useState(null);
    const [offset, setOffset] = useState({ x: 0, y: 0 });
    
    useEffect(() => {
      if (containerRef.current && nodes.length > 0) {
        const { width, height } = containerRef.current.getBoundingClientRect();
        const centerX = width / 2;
        const centerY = height / 2;
        const radius = Math.min(width, height) / 2 - 40;
        const newPositions = { ...positions };
        let needsUpdate = false;

        nodes.forEach((node, index) => {
          if (!newPositions[node]) {
            needsUpdate = true;
            const angle = (2 * Math.PI * index) / nodes.length;
            newPositions[node] = {
              x: centerX + radius * Math.cos(angle),
              y: centerY + radius * Math.sin(angle)
            };
          }
        });
        if(needsUpdate) {
          setPositions(newPositions);
        }
      }
    }, [nodes, containerRef.current?.clientWidth]);

    const handleMouseDown = (e, node) => {
      const pos = positions[node];
      const svgPoint = containerRef.current.createSVGPoint();
      svgPoint.x = e.clientX;
      svgPoint.y = e.clientY;
      const screenToSvg = e.target.getScreenCTM().inverse();
      const { x, y } = svgPoint.matrixTransform(screenToSvg);

      setDraggingNode(node);
      setOffset({ x: pos.x - x, y: pos.y - y });
    };

    const handleMouseMove = (e) => {
      if (!draggingNode) return;
      
      const svgPoint = containerRef.current.createSVGPoint();
      svgPoint.x = e.clientX;
      svgPoint.y = e.clientY;
      const screenToSvg = e.target.getScreenCTM().inverse();
      const { x, y } = svgPoint.matrixTransform(screenToSvg);
      
      setPositions(prev => ({
        ...prev,
        [draggingNode]: { x: x + offset.x, y: y + offset.y }
      }));
    };

    const handleMouseUp = () => {
      setDraggingNode(null);
    };

    const currentHighlightedNodes = currentStepData?.highlightIndices || [];

    return (
      <div className="bg-gradient-to-br from-gray-900 to-black rounded-xl p-6 min-h-[550px] border border-gray-700 shadow-inner">
        <div className="flex justify-between items-center mb-4">
          <h3 className="text-lg font-bold text-purple-400">🌐 Graph Canvas</h3>
          <div className="text-sm text-gray-400">
            Nodes: {nodes.length} | Edges: {edges.length}
          </div>
        </div>
        {nodes.length === 0 ? (
          <div className="text-gray-400 text-center py-10">
            <div className="text-4xl mb-4">🌐</div>
            <div className="text-lg mb-2">No graph data available.</div>
            <div className="text-sm">Add nodes and edges to begin.</div>
          </div>
        ) : (
          <svg width="100%" height="350px" ref={containerRef} onMouseMove={handleMouseMove} onMouseUp={handleMouseUp} onMouseLeave={handleMouseUp}>
            <defs>
              <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="9.5" refY="3.5" orient="auto">
                <polygon points="0 0, 10 3.5, 0 7" fill="#6B7280" />
              </marker>
            </defs>
            {edges.map((edge, index) => {
              const fromPos = positions[edge.from];
              const toPos = positions[edge.to];
              if (!fromPos || !toPos) return null;
              const isHighlighted = currentHighlightedNodes.includes(edge.from) && currentHighlightedNodes.includes(edge.to);
              return (
                <g key={`edge-${index}`}>
                  <line x1={fromPos.x} y1={fromPos.y} x2={toPos.x} y2={toPos.y} stroke={isHighlighted ? "#FBBF24" : "#4B5563"} strokeWidth="2" markerEnd={graphType === 'directed' ? "url(#arrowhead)" : ""} />
                  <text x={(fromPos.x + toPos.x) / 2} y={(fromPos.y + toPos.y) / 2} fill="#FBBF24" fontSize="12" textAnchor="middle" dy="-5">{edge.weight}</text>
                </g>
              );
            })}
            {nodes.map(node => {
              const pos = positions[node];
              if (!pos) return null;
              const isHighlighted = currentHighlightedNodes.includes(node);
              return (
                <g key={`node-${node}`} onMouseDown={(e) => handleMouseDown(e, node)} style={{ cursor: 'grab' }}>
                  <circle cx={pos.x} cy={pos.y} r="20" fill={isHighlighted ? "#8B5CF6" : "#3B82F6"} stroke={isHighlighted ? "#A78BFA" : "#1E40AF"} strokeWidth="3" className="transition-all" />
                  <text x={pos.x} y={pos.y} fill="white" fontSize="14" textAnchor="middle" dy="5" className="font-bold pointer-events-none">{node}</text>
                </g>
              );
            })}
          </svg>
        )}
      </div>
    );
  };

  // Navigate to previous step
  const goToPreviousStep = () => {
    if (currentStep > 0) {
      setCurrentStep(currentStep - 1);
      renderStep(currentStep - 1);
    }
  };

  // Navigate to next step
  const goToNextStep = () => {
    if (currentStep < traversalSteps.length - 1) {
      setCurrentStep(currentStep + 1);
      renderStep(currentStep + 1);
    }
  };

  // Auto play traversal steps
  const toggleAutoPlay = () => {
    if (autoPlayIntervalRef.current) {
      clearInterval(autoPlayIntervalRef.current);
      autoPlayIntervalRef.current = null;
    } else {
      autoPlayIntervalRef.current = setInterval(() => {
        setCurrentStep(prev => {
          const nextStep = prev + 1;
          if (nextStep < traversalSteps.length) {
            renderStep(nextStep);
            return nextStep;
          } else {
            clearInterval(autoPlayIntervalRef.current);
            autoPlayIntervalRef.current = null;
            return prev;
          }
        });
      }, 1000);
    }
  };

  // Render a specific step with enhanced highlighting
  const renderStep = (stepIndex) => {
    if (stepIndex < 0 || stepIndex >= traversalSteps.length) return;
    const step = traversalSteps[stepIndex];

    // Update array visualization
    if (step.arrayState) {
      setArray(step.arrayState);
    }

    // Reset all element states first
    resetElementStates();

    // Apply highlighting based on the step
    if (step.highlightIndices && step.highlightIndices.length > 0) {
      setTimeout(() => {
        const currentAlgorithmGroup = selectedAlgorithm.split(':')[0];
        if (currentAlgorithmGroup === 'sort' || currentAlgorithmGroup === 'search') {
          step.highlightIndices.forEach(index => {
            if (typeof index !== 'number') return;
            let state = 'active'; // Default highlight
            if (step.description.toLowerCase().includes('found')) state = 'found';
            else if (step.description.toLowerCase().includes('swap')) state = 'swap';
            else if (step.description.toLowerCase().includes('minimum')) state = 'compared';
            else if (step.description.toLowerCase().includes('sorted') || step.description.toLowerCase().includes('correct position')) state = 'ok';
            else if (step.description.toLowerCase().includes('pivot')) state = 'pivot';
            updateElementState(index, state);
          });
        }
      }, 50);
    }

    // Update stats
    setComparisons(step.comparisons || 0);
    setSwaps(step.swaps || 0);
    setIterations(step.iterations || 0);
  };

// Array Steps Panel component with scroll functionality
const ArrayStepsPanel = () => {
  // Auto-scroll to current step
  React.useEffect(() => {
    if (stepRefs.current[currentStep] && stepsContainerRef.current) {
      const currentStepElement = stepRefs.current[currentStep];
      const container = stepsContainerRef.current;
      
      // Calculate if the element is in view
      const containerRect = container.getBoundingClientRect();
      const elementRect = currentStepElement.getBoundingClientRect();
      
      // Check if element is outside the visible area
      if (elementRect.top < containerRect.top || elementRect.bottom > containerRect.bottom) {
        // Use custom slow scroll instead of native smooth scroll
        const targetScrollTop = currentStepElement.offsetTop - container.offsetTop - (container.clientHeight / 2) + (currentStepElement.clientHeight / 2);
        
        // Smooth scroll with custom duration
        const startScrollTop = container.scrollTop;
        const scrollDistance = targetScrollTop - startScrollTop;
        const duration = 800; // 800ms for slower scroll
        const startTime = performance.now();
        
        const animateScroll = (currentTime) => {
          const elapsed = currentTime - startTime;
          const progress = Math.min(elapsed / duration, 1);
          
          // Easing function for smoother animation
          const easeInOutQuad = progress => progress < 0.5 
            ? 2 * progress * progress 
            : 1 - Math.pow(-2 * progress + 2, 2) / 2;
          
          container.scrollTop = startScrollTop + (scrollDistance * easeInOutQuad(progress));
          
          if (progress < 1) {
            requestAnimationFrame(animateScroll);
          }
        };
        
        requestAnimationFrame(animateScroll);
      }
    }
  }, [currentStep]);

  return (
    <>
      <div className="flex justify-between items-center mb-4">
        <h3 className="text-xl font-bold text-green-400 flex items-center">
          📋 Algorithm Steps
        </h3>
        <div className="bg-gradient-to-r from-gray-700 to-gray-600 text-sm text-gray-200 px-3 py-1 rounded-full border border-gray-500">
          Step {currentStep + 1} of {traversalSteps.length}
        </div>
      </div>
      <div ref={stepsContainerRef} className="bg-gradient-to-b from-gray-900 to-black p-4 rounded-lg h-80 overflow-y-auto mb-4 border border-gray-700 shadow-inner scrollbar-thin scrollbar-track-gray-800 scrollbar-thumb-gray-600">
        {traversalSteps.length > 0 ? (
          traversalSteps.map((step, index) => (
            <div
              key={index}
              ref={el => stepRefs.current[index] = el}
              className={`p-4 rounded-lg mb-3 transition-all duration-300 cursor-pointer border-l-4 ${
                index === currentStep 
                  ? 'bg-gradient-to-r from-blue-900 to-blue-800 border-blue-400 shadow-lg transform scale-105' 
                  : 'bg-gradient-to-r from-gray-800 to-gray-750 border-gray-600 hover:from-gray-750 hover:to-gray-700 hover:border-gray-500'
              }`}
              onClick={() => {
                setCurrentStep(index);
                if (renderStep) renderStep(index);
              }}
            >
              <div className="flex justify-between items-start mb-2">
                <div className="font-semibold text-white">Step {index + 1}</div>
                {index === currentStep && <div className="text-blue-400 text-sm animate-pulse">👁️ Current</div>}
              </div>
              <div className="text-gray-100 mb-2 leading-relaxed">{step.description}</div>
              {step.arrayState && step.arrayState.length > 0 && (
                <div className="text-sm mt-2 p-2 bg-gray-900 rounded border border-gray-600">
                  <div className="font-bold text-yellow-400 mb-1">Current Array:</div>
                  <div className="font-mono text-green-400">
                    [{step.arrayState.join(', ')}]
                  </div>
                </div>
              )}
              <div className="text-xs mt-2 text-gray-400 flex space-x-4">
                <span>🔍 {step.comparisons} comparisons</span>
                <span>🔄 {step.swaps} swaps</span>
                <span>🔁 {step.iterations} iterations</span>
              </div>
            </div>
          ))
        ) : (
          <div className="text-gray-400 p-6 text-center">
            <div className="text-4xl mb-4">🚀</div>
            <div className="text-lg mb-2">Ready to visualize!</div>
            <div className="text-sm">Run an algorithm to see the steps.</div>
          </div>
        )}
      </div>
    </>
  );
};

  // Graph Steps Panel component
  const GraphStepsPanel = () => (
    <>
      <div className="flex justify-between items-center mb-4">
        <h3 className="text-xl font-bold text-green-400 flex items-center">
          📋 Graph Traversal Steps
        </h3>
        <div className="bg-gradient-to-r from-gray-700 to-gray-600 text-sm text-gray-200 px-3 py-1 rounded-full border border-gray-500">
          Step {currentStep + 1} of {traversalSteps.length}
        </div>
      </div>
      <div ref={stepsContainerRef} className="bg-gradient-to-b from-gray-900 to-black p-4 rounded-lg h-80 overflow-y-auto mb-4 border border-gray-700 shadow-inner scrollbar-thin scrollbar-track-gray-800 scrollbar-thumb-gray-600">
        {traversalSteps.length > 0 ? (
          traversalSteps.map((step, index) => (
            <div
              key={index}
              ref={el => stepRefs.current[index] = el}
              className={`p-4 rounded-lg mb-3 transition-all duration-300 cursor-pointer border-l-4 ${ index === currentStep ? 'bg-gradient-to-r from-purple-900 to-purple-800 border-purple-400 shadow-lg transform scale-105' : 'bg-gradient-to-r from-gray-800 to-gray-750 border-gray-600 hover:from-gray-750 hover:to-gray-700 hover:border-gray-500' }`}
              onClick={() => setCurrentStep(index)}
            >
              <div className="flex justify-between items-start mb-2">
                <div className="font-semibold text-white">Step {index + 1}</div>
                {index === currentStep && <div className="text-purple-400 text-sm animate-pulse">👁️ Current</div>}
              </div>
              <div className="text-gray-100 mb-2 leading-relaxed">{step.description}</div>
              {step.graphState?.distances && (
                <div className="text-sm mt-2 p-2 bg-gray-900 rounded border border-gray-600">
                  <div className={`font-bold ${step.graphState.final ? 'text-green-400' : 'text-yellow-400'}`}>
                    {step.graphState.final ? 'Final Distances:' : 'Current Distances:'}
                  </div>
                  <div className="font-mono text-gray-300 grid grid-cols-2 sm:grid-cols-3 gap-1 mt-1">
                    {Object.entries(step.graphState.distances).map(([node, dist]) => (
                      <div key={node}>{node}: <span className="text-white">{dist === Infinity ? '∞' : dist}</span></div>
                    ))}
                  </div>
                </div>
              )}
              <div className="text-xs mt-2 text-gray-400 flex space-x-4">
                <span>🔍 {step.comparisons} comparisons</span>
                <span>🔁 {step.iterations} iterations</span>
              </div>
              </div>
          ))
        ) : (
          <div className="text-gray-400 p-6 text-center">
            <div className="text-4xl mb-4">🌐</div>
            <div className="text-lg mb-2">Ready to visualize!</div>
            <div className="text-sm">Run a graph algorithm to see the steps.</div>
          </div>
        )}
      </div>
    </>
  );

  return (
    <div className="min-h-screen bg-gradient-to-b from-gray-900 to-gray-800 text-gray-100 pb-16">
      {/* Header */}
      <header className="sticky top-0 z-30 bg-gray-900 bg-opacity-90 backdrop-blur-md border-b border-gray-700 shadow-lg">
        <div className="container mx-auto px-4 py-4">
          <h1 className="text-3xl font-bold mb-3 bg-gradient-to-r from-blue-400 to-purple-500 bg-clip-text text-transparent">
            🚀 Enhanced Algorithm Visualizer
          </h1>
          <div className="flex flex-wrap gap-3 items-center mb-3">
            <select 
              value={selectedAlgorithm}
              onChange={(e) => setSelectedAlgorithm(e.target.value)}
              className="bg-gray-800 text-gray-100 border border-gray-600 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
              disabled={isRunning}
            >
              <optgroup label="🔄 Sorting Algorithms">
                <option value="sort:bubble">Bubble Sort</option>
                <option value="sort:selection">Selection Sort</option>
                <option value="sort:insertion">Insertion Sort</option>
                <option value="sort:merge">Merge Sort</option>
                <option value="sort:quick">Quick Sort</option>
                <option value="sort:heap">Heap Sort</option>
              </optgroup>
              <optgroup label="🔍 Search Algorithms">
                <option value="search:linear">Linear Search</option>
                <option value="search:binary">Binary Search</option>
                <option value="search:jump">Jump Search</option>
                <option value="search:interpolation">Interpolation Search</option>
              </optgroup>
              <optgroup label="🌳 Graph Algorithms">
                <option value="graph:dijkstra">Dijkstra's Algorithm</option>
                <option value="graph:bellmanford">Bellman-Ford Algorithm</option>
                <option value="graph:prims">Prim's Algorithm</option>
                <option value="graph:kruskal">Kruskal's Algorithm</option>
                <option value="graph:bfs">Graph BFS</option>
                <option value="graph:dfs">Graph DFS</option>
              </optgroup>
            </select>

            <div className="bg-gray-800 px-3 py-2 rounded-lg border border-gray-600">
              <label className="text-sm mr-2">⚡ Speed:</label>
              <input 
                type="range" 
                min="1" 
                max="200" 
                value={speed}
                onChange={(e) => setSpeed(parseInt(e.target.value))}
                className="w-20"
                disabled={isRunning}
              />
              <span className="text-xs ml-1">{speed}%</span>
            </div>

            <div className="bg-gray-800 px-3 py-2 rounded-lg border border-gray-600">
              <label className="text-sm mr-2">📊 Elements:</label>
              <input 
                type="number" 
                min="5" 
                max="20" 
                value={count}
                onChange={(e) => setCount(parseInt(e.target.value))}
                className="w-16 bg-gray-700 rounded px-1 text-center focus:ring-2 focus:ring-blue-500"
                disabled={isRunning}
              />
            </div>

            <div className="bg-gray-800 px-3 py-2 rounded-lg border border-gray-600">
              <label className="text-sm mr-2">🎯 Target:</label>
              <input 
                type="number" 
                value={target}
                onChange={(e) => setTarget(e.target.value)}
                placeholder="auto"
                className="w-20 bg-gray-700 rounded px-1 text-center focus:ring-2 focus:ring-blue-500"
                disabled={isRunning}
              />
            </div>

            <div className="bg-gray-800 px-3 py-2 rounded-lg border border-gray-600 w-[1000px]">
              <label className="text-sm mr-2">📝 Custom:</label>
              <input 
                type="text" 
                value={customArray}
                onChange={(e) => setCustomArray(e.target.value)}
                placeholder="5,3,8,1,9"
                className="w-200 bg-gray-700 rounded px-2 py-1 text-sm focus:ring-2 focus:ring-blue-500"
                disabled={isRunning}
              />
              <button 
                onClick={applyCustomArray}
                className="bg-blue-600 w-20 hover:bg-blue-300 hover:text-black px-2 py-1 rounded ml-2 text-sm transition-colors"
                disabled={isRunning}
              >
                Apply
              </button>
            </div>

            <div className="flex-grow"></div>

            <div className="flex gap-2">
              <button 
                onClick={generateArray}
                className="bg-gradient-to-r from-blue-600 to-blue-700 hover:from-blue-500 hover:to-blue-600 px-4 py-2 rounded-lg font-medium transition-all transform hover:scale-105 shadow-lg"
                disabled={isRunning}
              >
                🎲 Generate
              </button>

              <button 
                onClick={runAlgorithm}
                disabled={isRunning}
                className="bg-gradient-to-r from-green-600 to-green-700 hover:from-green-500 hover:to-green-600 px-4 py-2 rounded-lg font-medium disabled:opacity-50 transition-all transform hover:scale-105 shadow-lg"
              >
                {isRunning ? '🔄 Running...' : '▶️ Run'}
              </button>

              <button 
                onClick={togglePause}
                disabled={!isRunning}
                className="bg-gradient-to-r from-yellow-600 to-yellow-700 hover:from-yellow-500 hover:to-yellow-600 px-4 py-2 rounded-lg font-medium disabled:opacity-50 transition-all transform hover:scale-105 shadow-lg"
              >
                {isPaused ? '▶️ Resume' : '⏸️ Pause'}
              </button>

              <button 
                onClick={stopAlgorithm}
                className="bg-gradient-to-r from-red-600 to-red-700 hover:from-red-500 hover:to-red-600 px-4 py-2 rounded-lg font-medium transition-all transform hover:scale-105 shadow-lg"
              >
                ⏹️ Stop
              </button>
            </div>
          </div>

          <div className="flex items-center justify-between">
            <p className="text-sm text-gray-400">
              💡 <strong>Tip:</strong> Click and drag graph nodes to rearrange them!
            </p>
            
            {isRunning && (
              <div className="flex items-center space-x-2">
                <div className="animate-spin rounded-full h-4 w-4 border-2 border-green-400 border-t-transparent"></div>
                <span className="text-sm text-green-400 font-medium">Algorithm Running</span>
              </div>
            )}
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="container mx-auto px-4 py-6">
        <div className="flex flex-col lg:flex-row gap-6">
          {/* Visualization Panel */}
          <div className="lg:w-1/2 bg-gradient-to-br from-gray-800 to-gray-850 border border-gray-700 rounded-xl p-6 shadow-xl">
            {renderActiveVisualization()}
            {/* Enhanced Stats */}
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mt-6">
              <div className="bg-gradient-to-br from-green-800 to-green-900 p-4 rounded-lg border border-green-700 shadow-lg">
                <div className="text-3xl font-bold text-green-400 mb-1">{comparisons}</div>
                <div className="text-sm text-green-300">🔍 Comparisons</div>
              </div>
              <div className="bg-gradient-to-br from-blue-800 to-blue-900 p-4 rounded-lg border border-blue-700 shadow-lg">
                <div className="text-3xl font-bold text-blue-400 mb-1">{swaps}</div>
                <div className="text-sm text-blue-300">🔄 Swaps</div>
              </div>
              <div className="bg-gradient-to-br from-purple-800 to-purple-900 p-4 rounded-lg border border-purple-700 shadow-lg">
                <div className="text-3xl font-bold text-purple-400 mb-1">{iterations}</div>
                <div className="text-sm text-purple-300">🔁 Iterations</div>
              </div>
              <div className="bg-gradient-to-br from-orange-800 to-orange-900 p-4 rounded-lg border border-orange-700 shadow-lg">
                <div className="text-3xl font-bold text-orange-400 mb-1">{traversalSteps.length}</div>
                <div className="text-sm text-orange-300">📋 Steps</div>
              </div>
            </div>
          </div>

          {/* Enhanced Traversal Steps Panel */}
          <div className="lg:w-1/2 bg-gradient-to-br from-gray-800 to-gray-850 border border-gray-700 rounded-xl p-6 shadow-xl">
            {selectedAlgorithm.startsWith('graph:') ? <GraphStepsPanel/> : <ArrayStepsPanel/>}

            {/* Enhanced Control Buttons */}
            <div className="flex gap-2 flex-wrap mb-6">
              <button 
                onClick={goToPreviousStep}
                disabled={currentStep === 0 || traversalSteps.length === 0}
                className="bg-gradient-to-r from-gray-700 to-gray-800 hover:from-gray-600 hover:to-gray-700 px-4 py-2 rounded-lg disabled:opacity-50 transition-all transform hover:scale-105 shadow-lg"
              >
                ⏮️ Previous
              </button>
              <button 
                onClick={goToNextStep}
                disabled={currentStep === traversalSteps.length - 1 || traversalSteps.length === 0}
                className="bg-gradient-to-r from-gray-700 to-gray-800 hover:from-gray-600 hover:to-gray-700 px-4 py-2 rounded-lg disabled:opacity-50 transition-all transform hover:scale-105 shadow-lg"
              >
                Next ⏭️
              </button>
              <button 
                onClick={toggleAutoPlay}
                disabled={traversalSteps.length === 0}
                className="bg-gradient-to-r from-purple-700 to-purple-800 hover:from-purple-600 hover:to-purple-700 px-4 py-2 rounded-lg disabled:opacity-50 transition-all transform hover:scale-105 shadow-lg"
              >
                {autoPlayIntervalRef.current ? '⏸️ Pause Auto' : '▶️ Auto Play'}
              </button>
            </div>

            {/* Enhanced Code Panel */}
            <div className="bg-gradient-to-b from-gray-900 to-black rounded-lg border border-gray-700 shadow-inner">
              <div className="bg-gradient-to-r from-gray-800 to-gray-750 px-4 py-2 rounded-t-lg border-b border-gray-600">
                <h4 className="text-lg font-bold text-yellow-400 flex items-center">
                  💻 Algorithm Implementation
                </h4>
              </div>
              <pre className="p-4 overflow-x-auto text-sm max-h-64 text-green-400 font-mono leading-relaxed">
                {algorithmCode[selectedAlgorithm.split(':')[1]] || "// Select an algorithm to view its implementation\n// This shows the actual code structure"}
              </pre>
            </div>
          </div>
        </div>

        {/* Graph Controls and Visualization */}
        <div className="mt-6 bg-gradient-to-br from-gray-800 to-gray-850 border border-gray-700 rounded-xl p-6 shadow-xl">
          <h3 className="text-xl font-bold text-purple-400 mb-4">🌐 Graph Controls</h3>
          
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {/* Graph Controls */}
            <div className="space-y-4">
              <div className="bg-gray-900 p-4 rounded-lg border border-gray-700">
                <h4 className="text-lg font-bold text-blue-400 mb-3">Graph Controls</h4>
                
                <div className="grid grid-cols-2 gap-4 mb-4">
                  <div>
                    <label className="block text-sm text-gray-300 mb-1">Graph Type</label>
                    <select 
                      value={graphType}
                      onChange={(e) => setGraphType(e.target.value)}
                      className="w-full bg-gray-800 text-gray-100 border border-gray-600 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                    >
                      <option value="directed">Directed</option>
                      <option value="undirected">Undirected</option>
                    </select>
                  </div>
                  
                  <div>
                    <label className="block text-sm text-gray-300 mb-1">Start Node</label>
                    <select 
                      value={startNode}
                      onChange={(e) => setStartNode(e.target.value)}
                      className="w-full bg-gray-800 text-gray-100 border border-gray-600 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                    >
                      {graphNodes.map(node => (
                        <option key={node} value={node}>{node}</option>
                      ))}
                    </select>
                  </div>
                </div>
                
                <div className="mb-4">
                  <label className="block text-sm text-gray-300 mb-1">Graph JSON Input</label>
                  <textarea 
                    value={graphInput}
                    onChange={(e) => setGraphInput(e.target.value)}
                    placeholder='{"A": {"B": 1}, "B": {"A": 1}}'
                    className="w-full h-20 bg-gray-800 text-gray-100 border border-gray-600 rounded-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                  />
                  <button 
                    onClick={parseGraphInput}
                    className="mt-2 bg-blue-600 hover:bg-blue-500 px-3 py-1 rounded text-sm transition-colors w-full"
                  >
                    Parse Graph
                  </button>
                </div>
                
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <label className="block text-sm text-gray-300 mb-1">Add Node</label>
                    <div className="flex">
                      <input 
                        type="text" 
                        value={nodeValue}
                        onChange={(e) => setNodeValue(e.target.value)}
                        placeholder="Node name"
                        className="flex-1 bg-gray-800 text-gray-100 border border-gray-600 rounded-l-lg px-3 py-2 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      />
                      <button 
                        onClick={addNode}
                        className="bg-green-600 hover:bg-green-500 px-3 py-2 rounded-r-lg transition-colors"
                      >
                        +
                      </button>
                    </div>
                  </div>
                  
                  <div>
                    <label className="block text-sm text-gray-300 mb-1">Add Edge</label>
                    <div className="flex mb-1">
                      <input 
                        type="text" 
                        value={fromNode}
                        onChange={(e) => setFromNode(e.target.value)}
                        placeholder="From"
                        className="w-16 bg-gray-800 text-gray-100 border border-gray-600 rounded-l-lg px-2 py-1 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      />
                      <input 
                        type="text" 
                        value={toNode}
                        onChange={(e) => setToNode(e.target.value)}
                        placeholder="To"
                        className="w-16 bg-gray-800 text-gray-100 border border-gray-600 border-l-0 border-r-0 px-2 py-1 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      />
                      <input 
                        type="number" 
                        value={edgeWeight}
                        onChange={(e) => setEdgeWeight(parseInt(e.target.value) || 1)}
                        placeholder="Weight"
                        className="w-16 bg-gray-800 text-gray-100 border border-gray-600 rounded-r-lg px-2 py-1 focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                      />
                    </div>
                    <button 
                      onClick={addEdge}
                      className="bg-blue-600 hover:bg-blue-500 px-3 py-1 rounded text-sm transition-colors w-full"
                    >
                      Add Edge
                    </button>
                  </div>
                </div>
              </div>
              
              <div className="bg-gray-900 p-4 rounded-lg border border-gray-700">
                <h4 className="text-lg font-bold text-green-400 mb-3">Current Graph</h4>
                <div className="text-sm font-mono text-gray-300 max-h-40 overflow-y-auto">
                  {Object.keys(graph).length === 0 ? (
                    <div className="text-gray-500">No graph data</div>
                  ) : (
                    <pre>{JSON.stringify(graph, null, 2)}</pre>
                  )}
                </div>
              </div>
            </div>
            
            {/* Graph Visualization */}
            <div className="bg-gradient-to-b from-gray-900 to-black rounded-xl p-4 border border-gray-700 shadow-inner min-h-[400px]">
              <GraphCanvas 
                nodes={graphNodes}
                edges={graphEdges}
                positions={nodePositions}
                setPositions={setNodePositions}
                currentStepData={traversalSteps[currentStep]}
              />
            </div>
          </div>
        </div>
      </main>
    </div>
  );
};

export default AlgorithmVisualizer;