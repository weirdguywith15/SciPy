import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
from scipy.interpolate import griddata

# Set up the figure with multiple subplots
plt.style.use('dark_background')
fig = plt.figure(figsize=(18, 10))
fig.canvas.manager.set_window_title('Galaxy Energy-Phase Spiral with Space Tensor & Dimensional Flux')

# Create grid for subplots
gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], width_ratios=[2, 1])

# Create the main galaxy plot (2D spiral)
ax_galaxy = fig.add_subplot(gs[0, 0])
ax_galaxy.set_aspect('equal')
ax_galaxy.set_xlim(-1.1, 1.1)
ax_galaxy.set_ylim(-1.1, 1.1)
ax_galaxy.set_title("Galaxy Energy-Phase Spiral", color='white', fontsize=14)
ax_galaxy.axis('off')

# Create the 3D tensor heatmap
ax_tensor = fig.add_subplot(gs[0, 1], projection='3d')
ax_tensor.set_title("Space Tensor & Dimensional Flux", color='white', fontsize=14)
ax_tensor.set_facecolor('black')
ax_tensor.xaxis.pane.fill = False
ax_tensor.yaxis.pane.fill = False
ax_tensor.zaxis.pane.fill = False
ax_tensor.grid(False)

# Set initial parameters
params = {
    'rotation_speed': 1.0,
    'phase_shift': 0.5,
    'energy_level': 1.0,
    'anomaly_strength': 0.3,
    'arm_count': 4
}

# Define gravitational anomalies
anomalies = [
    {'x': 0.3, 'y': 0.7, 'strength': 0.8, 'size': 0.15},
    {'x': -0.5, 'y': -0.2, 'strength': 0.6, 'size': 0.1},
    {'x': 0.6, 'y': -0.4, 'strength': 0.7, 'size': 0.2}
]

# Text area for anomaly status
anomaly_text = fig.text(0.75, 0.85, 'Anomaly Status: Stable', fontsize=12, color='#4CAF50',
                        bbox=dict(facecolor='black', alpha=0.7, edgecolor='gray', boxstyle='round,pad=0.5'))

# Create sliders
ax_sliders = []
sliders = {}

slider_props = [
    ('rotation_speed', 'Rotation Speed', 0.1, 2.0, params['rotation_speed']),
    ('phase_shift', 'Phase Shift', 0.0, 1.0, params['phase_shift']),
    ('energy_level', 'Energy Level', 0.5, 2.0, params['energy_level']),
    ('anomaly_strength', 'Anomaly Strength', 0.0, 1.0, params['anomaly_strength']),
    ('arm_count', 'Spiral Arms', 2, 8, params['arm_count'])
]

# Add sliders for each parameter
for i, (param_name, label, min_val, max_val, init_val) in enumerate(slider_props):
    ax = fig.add_subplot(gs[1, 0], position=[0.1, 0.1 - i*0.05, 0.65, 0.03])
    ax_sliders.append(ax)
    step = 0.1 if param_name != 'arm_count' else 1
    sliders[param_name] = Slider(
        ax=ax, 
        label=label, 
        valmin=min_val, 
        valmax=max_val, 
        valinit=init_val,
        valstep=step,
        color='#3498db'
    )

# Get color for energy level
def get_energy_color(energy):
    """Map energy level to a color from blue (low) to yellow (high)"""
    energy = max(0, min(1, energy)) # Normalize between 0 and 1
    
    cmap = cm.get_cmap('plasma')
    return cmap(energy)

# Calculate gravitational anomaly effect
def calculate_anomaly_effect(x, y):
    """Calculate the effect of gravitational anomalies at a given point"""
    effect = {'distortion': 0, 'energy_shift': 0, 'tensor_field': 0}
    
    for anomaly in anomalies:
        dx = x - anomaly['x']
        dy = y - anomaly['y']
        distance = np.sqrt(dx**2 + dy**2)
        
        if distance < anomaly['size'] * 2:
            strength = (1 - distance / (anomaly['size'] * 2)) * anomaly['strength'] * params['anomaly_strength']
            effect['distortion'] += strength * 0.5
            effect['energy_shift'] += strength * 0.6
            effect['tensor_field'] += strength * 0.8
    
    return effect

# Update anomaly status display
def update_anomaly_status():
    """Update the anomaly status display based on current anomaly strength"""
    strength = params['anomaly_strength']
    
    if strength < 0.2:
        status = "Stable"
        color = '#4CAF50'
    elif strength < 0.5:
        status = "Minor Anomalies"
        color = '#FFC107'
    elif strength < 0.8:
        status = "Significant Distortion"
        color = '#FF9800'
    else:
        status = "Critical Instability"
        color = '#F44336'
    
    anomaly_text.set_text(f'Anomaly Status: {status}')
    anomaly_text.set_color(color)

# Generate data for tensor field grid
def generate_tensor_field(time):
    """Generate the tensor field and dimensional flux data"""
    grid_size = 15
    x = np.linspace(-1, 1, grid_size)
    y = np.linspace(-1, 1, grid_size)
    z = np.linspace(-1, 1, grid_size)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    tensor_values = np.zeros((grid_size, grid_size, grid_size))
    
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                # Distance from center
                r = np.sqrt(X[i,j,k]2 + Y[i,j,k]**2 + Z[i,j,k]**2)
                
                if r > 0.9:  # Skip points outside our region of interest
                    continue
                
                # Base tensor value
                base_value = np.exp(-r*3) * (1 + 0.2*np.sin(r*10 + time*2))
                
                # Phase component
                phase = params['phase_shift'] * np.sin(r*6 + time)
                
                # Anomaly influence in 3D
                anomaly_effect = calculate_anomaly_effect(X[i,j,k], Y[i,j,k])
                tensor_field = base_value + anomaly_effect['tensor_field'] + phase
                
                # Energy level influence
                tensor_values[i,j,k] = tensor_field * params['energy_level']
    
    return X, Y, Z, tensor_values

# Initialize plot elements
scatter_galaxy = None
tensor_plot = None
time = 0

# Update visualization
def update_plot(frame):
    global scatter_galaxy, tensor_plot, time
    
    # Update time based on rotation speed
    time += 0.05 * params['rotation_speed']
    
    # Clear previous plots
    ax_galaxy.clear()
    ax_tensor.clear()
    
    # Set up galaxy plot
    ax_galaxy.set_aspect('equal')
    ax_galaxy.set_xlim(-1.1, 1.1)
    ax_galaxy.set_ylim(-1.1, 1.1)
    ax_galaxy.set_title("Galaxy Energy-Phase Spiral", color='white', fontsize=14)
    ax_galaxy.axis('off')
    
    # Set up tensor plot
    ax_tensor.set_title("Space Tensor & Dimensional Flux", color='white', fontsize=14)
    ax_tensor.set_facecolor('black')
    ax_tensor.set_xlim(-1, 1)
    ax_tensor.set_ylim(-1, 1)
    ax_tensor.set_zlim(-1, 1)
    ax_tensor.set_xlabel("X", color='white')
    ax_tensor.set_ylabel("Y", color='white')
    ax_tensor.set_zlabel("Z", color='white')
    ax_tensor.grid(False)
    
    # Generate particles for galaxy
    arm_count = int(params['arm_count'])
    num_particles = 5000
    arms = np.random.randint(0, arm_count, num_particles)
    
    # Random distribution biased toward outer regions
    r = np.sqrt(np.random.random(num_particles))
    arm_tightness = 2.5 * (1 - params['phase_shift'] * 0.5)
    
    # Rotation speed varies with radius
    rotation_speed = 1 - (r * 0.5)
    time_offset = time * rotation_speed
    
    # Calculate spiral arm positions
    theta = np.random.random(num_particles) * 2 * np.pi
    arm_offset = (arms / arm_count) * 2 * np.pi
    spiral_offset = r * arm_tightness * 2 * np.pi
    final_angle = theta + arm_offset + spiral_offset + time_offset
    
    # Calculate base positions
    x = np.cos(final_angle) * r
    y = np.sin(final_angle) * r
    
    # Apply phase shift
    phase_angle = params['phase_shift'] * np.sin(r * np.pi * 4 + time)
    phase_x = np.cos(final_angle + phase_angle) * (r * 0.1)
    phase_y = np.sin(final_angle + phase_angle) * (r * 0.1)
    
    x += phase_x
    y += phase_y
    
    # Apply gravitational anomalies
    energy_values = np.z
eros(num_particles)
    
    for i in range(num_particles):
        # Calculate anomaly effects
        anomaly_effect = calculate_anomaly_effect(x[i], y[i])
        
        # Apply distortion
        distortion_angle = np.arctan2(y[i], x[i]) + anomaly_effect['distortion'] * np.pi
        distortion_magnitude = anomaly_effect['distortion'] * r[i] * 0.3
        x[i] += np.cos(distortion_angle) * distortion_magnitude
        y[i] += np.sin(distortion_angle) * distortion_magnitude
        
        # Calculate energy
        base_energy = params['energy_level'] * (1 - r[i] * 0.5)
        energy_variation = np.sin(r[i] * 20 + arms[i] + time * 2) * 0.15
        energy_values[i] = base_energy + energy_variation + anomaly_effect['energy_shift']
    
    # Normalize energy values
    energy_values = np.clip(energy_values, 0, 1)
    
    # Set size based on energy
    sizes = 5 + energy_values * 15
    
    # Get colors based on energy
    colors_array = np.array([get_energy_color(e) for e in energy_values])
    
    # Plot the particles
    scatter_galaxy = ax_galaxy.scatter(x, y, s=sizes, c=colors_array, alpha=0.7)
    
    # Draw galaxy core
    core = plt.Circle((0, 0), 0.15, color='yellow', alpha=0.7)
    ax_galaxy.add_artist(core)
    
    # Draw anomalies if they're strong enough
    if params['anomaly_strength'] > 0.15:
        for anomaly in anomalies:
            visible_strength = params['anomaly_strength'] * anomaly['strength']
            if visible_strength > 0.2:
                anomaly_circle = plt.Circle(
                    (anomaly['x'], anomaly['y']), 
                    anomaly['size'], 
                    color='purple', 
                    alpha=visible_strength * 0.4)
                ax_galaxy.add_artist(anomaly_circle)
    
    # Generate and plot tensor field
    X, Y, Z, tensor_values = generate_tensor_field(time)
    
    # Create mask for non-zero values
    mask = tensor_values > 0.05
    
    # Create Isosurface at multiple levels
    levels = np.linspace(0.1, 0.8, 4)
    
    for level in levels:
        verts, faces, _, _ = measure_marching_cubes(tensor_values, level)
        
        # Scale to match our coordinate system
        verts = verts / (len(tensor_values) - 1) * 2 - 1
        
        # Plot the isosurface
        mesh = ax_tensor.plot_trisurf(
            verts[:, 0], verts[:, 1], verts[:, 2],
            triangles=faces,
            cmap='plasma',
            alpha=0.3,
            edgecolor='none'
        )
    
    # Plot some points to show the dimensional flux
    points_count = 100
    # Random points in a sphere
    phi = np.random.uniform(0, 2*np.pi, points_count)
    costheta = np.random.uniform(-1, 1, points_count)
    theta = np.arccos(costheta)
    r = np.random.uniform(0, 0.9, points_count)(1/3)
    
    x_points = r * np.sin(theta) * np.cos(phi)
    y_points = r * np.sin(theta) * np.sin(phi)
    z_points = r * np.cos(theta)
    
    # Calculate flux values at each point
    flux_values = np.zeros(points_count)
    for i in range(points_count):
        # Interpolate tensor value at this point
        base_r = np.sqrt(x_points[i]**2 + y_points[i]**2 + z_points[i]**2)
        base_flux = np.exp(-base_r*3) * (1 + 0.3*np.sin(base_r*8 + time*3))
        
        # Add anomaly effects (projected to this 3D point)
        xy_dist = np.sqrt(x_points[i]**2 + y_points[i]**2)
        if xy_dist > 0:
            scaled_x = x_points[i] / xy_dist
            scaled_y = y_points[i] / xy_dist
            anomaly_effect = calculate_anomaly_effect(scaled_x, scaled_y)
            flux_values[i] = base_flux + anomaly_effect['tensor_field'] * params['energy_level']
        else:
            flux_values[i] = base_flux
    
    # Normalize flux values
    flux_values = (flux_values - flux_values.min()) / (flux_values.max() - flux_values.min() + 1e-6)
    
    # Plot the flux points
    ax_tensor.scatter(
        x_points, y_points, z_points, 
        c=flux_values, 
        cmap='plasma', 
        s=10 + 20*flux_values, 
        alpha=0.6
    )
    
    # Update anomaly status
update_anomaly_status()
    
    return scatter_galaxy,

# Function to handle slider changes
def update_from_slider(val):
    for param_name, slider in sliders.items():
        if param_name == 'arm_count':
            params[param_name] = int(slider.val)
        else:
            params[param_name] = slider.val

# Connect sliders to update function
for param_name, slider in sliders.items():
    slider.on_changed(update_from_slider)

# Helper function for marching cubes algorithm (simplified version)
def measure_marching_cubes(volume, level):
    """A simplified implementation of marching cubes for isosurface extraction"""
    # Get the indices where values are above the level
    xx, yy, zz = np.where(volume > level)
    
    # Create vertices (we'll simplify by using cube corners)
    verts = np.array([[x-0.5, y-0.5, z-0.5] for x, y, z in zip(xx, yy, zz)])
    
    # Create simple faces (we'll make small tetrahedrons)
    faces = []
    for i in range(len(verts) - 3):
        if i % 4 == 0:  # Create a tetrahedron every 4 points
            faces.append([i, i+1, i+2])
            faces.append([i, i+2, i+3])
            faces.append([i, i+1, i+3])
            faces.append([i+1, i+2, i+3])
    
    # If not enough points, create a single tetrahedral
    if len(faces) == 0 and len(verts) >= 4:
        faces = [[0, 1, 2], [0, 2, 3], [0, 1, 3], [1, 2, 3]]
    
    return verts, np.array(faces), None, None

# Create the animation
ani = animation.FuncAnimation(fig, update_plot, frames=100, interval=50, blit=False)

# Add info text
info_text = """
Galaxy Energy-Phase Diagram with Space Tensor & Dimensional Flux

This visualization shows the relationship between energy distribution 
and phase shifting in a galactic structure, with gravitational anomalies 
creating distortions in spacetime.

The 3D plot displays the space tensor field and dimensional flux metrics,
revealing how gravitational anomalies affect the underlying structure of 
spacetime around the galaxy.

Adjust the sliders to explore different parameter configurations.
"""

fig.text(0.7, 0.35, info_text, fontsize=10, color='white',
         bbox=dict(facecolor='black', alpha=0.7, edgecolor='gray', boxstyle='round,pad=0.5'))

plt.tight_layout()
plt.show()
