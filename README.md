## Security Officer Shift Allocation

### Project Overview

The Security Officer Shift Allocation project aims to efficiently assign security officers to various companies based on their shift preferences and the companies' requirements using an efficient algorithm. This ensures optimal resource utilization and adherence to work constraints.

### Key Features

- **Preference-Based Allocation:** Considers security officers' preferred shifts to maximize satisfaction.
- **Constraint Handling:** Ensures officers work within specified minimum and maximum shift limits.
- **Daily Consistency:** Maintains the required number of officers for each company shift every day.
- **Algorithmic Approach:** Utilizes the Ford-Fulkerson method to solve the allocation problem effectively.

### Usefulness

This project demonstrates an ability to solve complex scheduling problems using algorithmic solutions. It is particularly useful for companies needing to manage workforce allocation efficiently, ensuring that staffing requirements are met while respecting employees' preferences and legal constraints.

### Basic Information

- **Inputs:** Preferences of security officers, company shift requirements, minimum and maximum shifts per officer.
- **Outputs:** An allocation plan satisfying all constraints or an indication that no valid allocation exists.
- **Complexity:** Optimized for practical use with a focus on efficient computation and resource management.

### Function Definition

```python
def allocate(preferences, officers_per_org, min_shifts, max_shifts):
    """
    Allocates security officers to companies based on preferences and requirements.

    Parameters:
    - preferences: List of lists. preferences[i][k] is 1 if officer SOi prefers shift Sk, else 0.
    - officers_per_org: List of lists. officers_per_org[j][k] is the number of officers company Cj needs for shift Sk.
    - min_shifts: Minimum number of shifts an officer must work in the month.
    - max_shifts: Maximum number of shifts an officer can work in the month.

    Returns:
    - A list of lists (allocation) if a valid assignment exists, otherwise `None`.
    """
```

### Input Description

1. **preferences:** `[[1, 1, 1], [1, 0, 1], [1, 1, 0], [0, 1, 0]]`
   - Binary values indicating officer shift preferences.

2. **officers_per_org:** `[[1, 0, 2], [0, 1, 0]]`
   - Number of officers required by each company for each shift.

3. **min_shifts:** `1`
   - Minimum number of shifts per officer.

4. **max_shifts:** `5`
   - Maximum number of shifts per officer.

### Output

- **Valid Allocation:** A nested list structure where `allocation[i][j][d][k]` is 1 if officer SOi is allocated to company Cj for shift Sk on day Dd.
- **No Valid Allocation:** Returns `None`.


### Running the File

To run the allocation function, use the following example code:

```python
# Example inputs
preferences = [[1, 1, 1], [1, 0, 1], [1, 1, 0], [0, 1, 0]]
officers_per_org = [[1, 0, 2], [0, 1, 0]]
min_shifts = 1
max_shifts = 5

# Call the allocate function
allocation = allocate(preferences, officers_per_org, min_shifts, max_shifts)

# Output the result
print(allocation)
```

This example demonstrates how to call the `allocate` function with sample inputs to generate and print the allocation plan.

### Complexity Analysis

- **Time Complexity:** O(m * n * n), where m is the number of companies and n is the number of security officers.
- **Space Complexity:** O(m * n), considering the auxiliary space for storing edges and vertices in the flow network

