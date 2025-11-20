"""
EnhancedSimulation class extending Simulation with enhanced DAG workflow support.

Provides error handlers, conditionals, machine filtering, and loops
for DAG workflows while maintaining backward compatibility.
"""

from typing import Dict, Set, List, Optional, Callable, Any, Tuple
import os
from datetime import datetime

from simulation import Simulation
from error_handler import ErrorHandler


class EnhancedSimulation(Simulation):
    """
    Extended Simulation class with enhanced DAG workflow support.
    
    Inheritance Pattern:
    - Inherits from Simulation using standard Python inheritance
    - Calls super().__init__(**kwargs) after extracting enhanced-specific kwargs
    - Adds enhanced attributes: attempt_number, iteration_count, error_handlers,
      condition, required_machine, allowed_machines, loop_enabled, etc.
    
    Method Overrides:
    
    1. __init__():
       Base: Initializes standard simulation attributes (path, job, dependencies, etc.)
       Override: Extracts enhanced kwargs first, initializes enhanced state, calls super(),
                 then performs machine filtering and definition-time conditional checks.
                 Early returns if marked fake (machine incompatible or conditional failed).
    
    2. set_directories():
       Base: Sets locdir, remdir, resdir based on path
       Override: Generates attempt-specific identifier/path for resubmissions before
                 calling super().set_directories()
    
    3. save_attempt():
       Base: Saves attempt information (if implemented in base)
       Override: Calls super() then stores enhanced metadata including timestamp,
                 attempt_number, iteration_count, error messages in attempt_history
    
    4. resubmit():
       Base: Resubmits simulation (if implemented in base)
       Override: Saves current attempt, increments attempt_number, generates new
                 path/identifier, resets indicators, updates directories/files.
                 Preserves dependencies automatically (stored by simid).
    
    5. check_status():
       Base: Checks simulation status and sets failed flag
       Override: Calls super() then invokes error handlers if failed. Can
                 trigger automatic resubmission with error handling by resubmit().
    
    6. progress():
       Base: Progresses simulation through stages (setup, submit, analyze, etc.)
             and triggers dependent simulations when finished
       Override: Checks execution-time conditionals first (skips if false), manages
                 loop conditions/limits, then calls super().progress(). After parent
                 call, triggers next_iteration() if loop enabled and finished.
    
    New Methods (Not in Base):
    - is_machine_compatible(): Checks if simulation can run on current machine
    - evaluate_condition(): Evaluates conditional function (returns bool)
    - skip(): Marks simulation as skipped (sets skipped=True, block=True)
    - generate_resubmission_path(): Generates new path/identifier for resubmission
    - next_iteration(): Prepares simulation for next loop iteration (increments count, resets indicators)
    - depends(): Override with strategy and selector parameters for merge strategies ('all', 'first')
    
    New Attributes (Not in Base):
    - attempt_number: Tracks resubmission attempts (0 = first run)
    - base_identifier/base_path: Original identifier/path before resubmissions
    - iteration_count: Current loop iteration (0 = first)
    - error_handlers: List of error handler objects/functions
    - attempt_history: List of attempt metadata dictionaries
    - condition: Conditional function (Callable[[Simulation], bool])
    - conditional_enabled: Whether conditional execution is active
    - skipped: Whether simulation was skipped
    - required_machine/allowed_machines: Machine filtering constraints
    - loop_enabled: Whether loop support is active
    - loop_condition: Function that determines if loop should continue
    - loop_max_iterations: Maximum loop iterations allowed
    - loop_modifier: Function called before each iteration to modify input/state
    - loop_variables: Dict for storing loop state
    - merge_strategy: Merge strategy for multiple dependencies ('all', 'first')
    - merge_selector: Selector function for merge strategy
    - _completed_dependencies: Set tracking completed dependencies for merge strategies
    
    Extension Strategy:
    - Pre-processing: Extracts enhanced kwargs, validates machine compatibility, evaluates
                      definition-time conditionals before calling super().__init__()
    - Post-processing: Adds error handler invocation after check_status(), loop iteration
                       handling after progress()
    - Side-by-side: Adds new functionality (machine filtering, loops, conditionals)
                    without modifying base behavior
    - Selective override: Only overrides methods needed for enhanced features
    
    Backward Compatibility:
    - DAG simulations unchanged: If no enhanced features used (no error_handlers, condition,
      required_machine, loop_enabled), behavior matches base Simulation class
    - Base class methods called via super() to preserve original functionality
    - New features opt-in: Enhanced features only active when explicitly configured
    - make_enhanced() function: Converts any Simulation to EnhancedSimulation dynamically,
      preserving all original attributes and behavior
    """
    
    # Allowed inputs extended from Simulation
    enhanced_allowed_inputs = set([
        'error_handlers', 'condition', 
        'required_machine', 'allowed_machines', 'loop_enabled',
        'loop_condition', 'loop_max_iterations', 'loop_modifier',
        'merge_strategy', 'merge_selector'
    ])
    
    def __init__(self, **kwargs):
        # Extract enhanced-specific kwargs
        enhanced_kwargs = {}
        for key in list(kwargs.keys()):
            if key in self.enhanced_allowed_inputs:
                enhanced_kwargs[key] = kwargs.pop(key)
        
        # Initialize enhanced state
        self.attempt_number: int = 0
        self.base_identifier: Optional[str] = None
        self.base_path: Optional[str] = None
        self.iteration_count: int = 0
        self.error_handlers: List[ErrorHandler] = list(enhanced_kwargs.get('error_handlers', []))
        self.attempt_history: List[Dict[str, Any]] = []
        self.condition: Optional[Callable[[Any], bool]] = enhanced_kwargs.get('condition', None)
        self.conditional_enabled: bool = self.condition is not None
        self.skipped: bool = False
        self.required_machine: Optional[str] = enhanced_kwargs.get('required_machine', None)
        self.allowed_machines: Set[str] = set(enhanced_kwargs.get('allowed_machines', []))
        self.loop_enabled: bool = enhanced_kwargs.get('loop_enabled', False)
        self.loop_condition: Optional[Callable[[Any], bool]] = enhanced_kwargs.get('loop_condition', None)
        self.loop_max_iterations: Optional[int] = enhanced_kwargs.get('loop_max_iterations', None)
        self.loop_modifier: Optional[Callable[[Any], None]] = enhanced_kwargs.get('loop_modifier', None)
        self.loop_variables: Dict[str, Any] = {}
        self.merge_strategy: str = enhanced_kwargs.get('merge_strategy', 'all')  # 'all', 'first'
        self.merge_selector: Optional[Callable[[List[Simulation]], Simulation]] = enhanced_kwargs.get('merge_selector', None)
        self._completed_dependencies: Set[int] = set()  # Track completed dependencies for merge strategies
        self._selected_dependency_id: Optional[int] = None  # Store selected dependency ID when merge_selector is used
        
        # Check for decorator-based configuration
        cls = self.__class__
        if hasattr(cls, '_required_machine'):
            self.required_machine = cls._required_machine
        if hasattr(cls, '_allowed_machines'):
            self.allowed_machines = cls._allowed_machines
        if hasattr(cls, '_loop_enabled'):
            self.loop_enabled = cls._loop_enabled
        if hasattr(cls, '_loop_condition'):
            self.loop_condition = cls._loop_condition
        if hasattr(cls, '_loop_max_iterations'):
            self.loop_max_iterations = cls._loop_max_iterations
        if hasattr(cls, '_loop_modifier'):
            self.loop_modifier = cls._loop_modifier
        if hasattr(cls, '_conditional_func'):
            self.condition = cls._conditional_func
            self.conditional_enabled = True
        
        # Call parent __init__ first to initialize base attributes
        super().__init__(**kwargs)
        
        # Machine filtering at definition time (after parent init so job is available)
        if self.required_machine or self.allowed_machines:
            if not self.is_machine_compatible():
                self.fake_sim = True  # Mark as fake to exclude from workflow
                return
        
        # Definition-time conditional check
        if self.conditional_enabled and hasattr(cls, '_conditional_execution_time'):
            if not cls._conditional_execution_time:  # Definition-time
                if not self.evaluate_condition():
                    self.fake_sim = True
                    return
        
        # Store base identifier and path for resubmissions
        if self.base_identifier is None:
            self.base_identifier = self.identifier
        if self.base_path is None:
            self.base_path = self.path
    
    def is_machine_compatible(self) -> bool:
        """Check if simulation is compatible with current machine."""
        # If required_machine or allowed_machines is set, we must check compatibility
        required = getattr(self, 'required_machine', None)
        allowed = getattr(self, 'allowed_machines', None)
        
        # If no machine restrictions, always compatible
        if not required and (not allowed or len(allowed) == 0):
            return True
        
        # Get current machine name from multiple sources
        machine_name = None
        
        # First try: job instance machine
        if hasattr(self, 'job') and self.job is not None:
            machine_name = getattr(self.job, 'machine', None)
        
        # Second try: Job class variable
        if machine_name is None:
            from machines import Job
            machine_name = Job.machine
        
        # Third try: ProjectManager machine
        if machine_name is None:
            from project_manager import ProjectManager
            if ProjectManager.machine is not None:
                machine_name = getattr(ProjectManager.machine, 'name', None)
        
        # If still no machine, can't verify - return False to be safe
        if machine_name is None:
            return False
        
        # Normalize machine names for comparison (strip whitespace, lowercase)
        machine_name = str(machine_name).strip().lower() if machine_name else None
        required = str(required).strip().lower() if required else None
        allowed = {str(a).strip().lower() for a in allowed} if allowed else set()
        
        # Check required_machine if set
        if required:
            compatible = machine_name == required
            if not compatible:
                # Explicitly mark as fake if incompatible
                self.fake_sim = True
            return compatible
        
        # Check allowed_machines if set
        if allowed and len(allowed) > 0:
            compatible = machine_name in allowed
            if not compatible:
                # Explicitly mark as fake if incompatible
                self.fake_sim = True
            return compatible
        
        return True
    
    def evaluate_condition(self) -> bool:
        """Evaluate conditional function."""
        if self.condition is None:
            return True
        try:
            return self.condition(self)
        except Exception as e:
            self.warn(f'Error evaluating condition: {e}')
            return False
    
    def get_conditional_info(self) -> str:
        """
        Get a human-readable description of the conditional configuration.
        
        Returns:
            String describing the conditional status and type
        """
        if not self.conditional_enabled:
            return ""
        
        # Handle edge case: conditional_enabled but condition is None
        if self.condition is None:
            return "cond:None[invalid]"
        
        # Check if we have execution-time vs definition-time info
        exec_time = "dynamic" if hasattr(self.__class__, '_conditional_execution_time') and self.__class__._conditional_execution_time else "static"
        
        # Get conditional description (with error handling)
        try:
            from conditionals import get_conditional_description
            cond_desc = get_conditional_description(self.condition)
            
            # Truncate very long descriptions to keep status line readable
            max_desc_len = 40
            if len(cond_desc) > max_desc_len:
                cond_desc = cond_desc[:max_desc_len-3] + "..."
        except Exception:
            cond_desc = "custom"
        
        # Check if skipped
        if self.skipped:
            return f"cond:{cond_desc}[{exec_time},skipped]"
        
        # Return conditional info without evaluation status (evaluation can be expensive/risky)
        return f"cond:{cond_desc}[{exec_time}]"
    
    def skip(self):
        """Mark simulation as skipped."""
        self.skipped = True
        self.block = True
    
    def set_directories(self):
        """Override to support resubmissions with new identifiers."""
        # Generate identifier with attempt number if resubmission
        if self.attempt_number > 0:
            self.identifier = f"{self.base_identifier}_attempt{self.attempt_number}"
            if self.base_path:
                self.path = f"{self.base_path}/attempt_{self.attempt_number}"
        
        # Call parent method
        super().set_directories()
    
    def generate_resubmission_path(self) -> Tuple[str, str]:
        """
        Generate new path and identifier for resubmission.
        
        Returns:
            Tuple of (new_path, new_identifier)
        """
        new_attempt = self.attempt_number + 1
        new_identifier = f"{self.base_identifier}_attempt{new_attempt}"
        new_path = f"{self.base_path}/attempt_{new_attempt}" if self.base_path else f"attempt_{new_attempt}"
        return new_path, new_identifier
    
    def save_attempt(self):
        """Enhanced save_attempt that stores more metadata."""
        # Call parent method
        super().save_attempt()
        
        # Store additional metadata
        attempt_data = {
            'timestamp': datetime.now().isoformat(),
            'attempt_number': self.attempt_number,
            'identifier': self.identifier,
            'path': self.path,
            'failed': self.failed,
            'finished': self.finished,
            'iteration_count': self.iteration_count,
        }
        
        # Try to capture error information
        if self.failed:
            attempt_data['error_type'] = 'simulation_failed'
            if hasattr(self, 'errfile') and self.errfile:
                errfile_path = os.path.join(self.locdir, self.errfile)
                if os.path.exists(errfile_path):
                    try:
                        with open(errfile_path, 'r', encoding='utf-8', errors='replace') as f:
                            attempt_data['error_message'] = f.read()[:1000]  # Limit size
                    except:
                        pass
        
        self.attempt_history.append(attempt_data)
    
    def resubmit(self):
        """
        Resubmit simulation after error.
        
        Increments attempt_number, generates new identifier/path,
        resets indicators, and preserves dependencies.
        """
        # Save current attempt
        self.save_attempt()
        
        # Increment attempt number
        self.attempt_number += 1
        
        # Generate new path and identifier
        new_path, new_identifier = self.generate_resubmission_path()
        self.path = new_path
        self.identifier = new_identifier
        
        # Reset indicators
        self.reset_indicators()
        
        # Update directories and files
        self.set_directories()
        self.set_files()
        
        # Dependencies are preserved automatically (they're stored by simid)
    
    def check_status(self):
        """Override to call error handlers on failure."""
        # Call parent method
        super().check_status()
        
        # If failed, call error handlers
        if self.failed and self.error_handlers:
            error_info = {
                'error_type': 'simulation_failed',
                'attempt_number': self.attempt_number,
                'iteration_count': self.iteration_count,
                'identifier': self.identifier,
            }
            
            should_resubmit = False
            handler_results = []  # Track which handlers were called and their results
            
            for handler in self.error_handlers:
                try:
                    # Check compatibility if handler supports it
                    if hasattr(handler, 'is_compatible'):
                        if not handler.is_compatible(self):
                            handler_name = getattr(handler, '__class__', type(handler)).__name__
                            self.warn(f'Error handler {handler_name} is not compatible with {self.__class__.__name__}, skipping')
                            handler_results.append((handler_name, 'skipped', 'incompatible'))
                            continue
                    
                    # Call the handler
                    if hasattr(handler, 'handle_error'):
                        result = handler.handle_error(self, error_info)
                    elif callable(handler):
                        result = handler(self, error_info)
                    else:
                        handler_results.append((str(handler), 'skipped', 'not_callable'))
                        continue
                    
                    handler_name = getattr(handler, '__class__', type(handler)).__name__
                    handler_results.append((handler_name, 'retry' if result else 'no_retry', result))
                    
                    # First handler that recommends retry triggers resubmission
                    # This implements a "first handler wins" strategy
                    # Handlers are tried in order, and the first one that says "retry" wins
                    if result:
                        should_resubmit = True
                        self.log(f'Error handler {handler_name} recommended retry, resubmitting', n=2)
                        break
                except Exception as e:
                    handler_name = getattr(handler, '__class__', type(handler)).__name__
                    self.warn(f'Error in error handler {handler_name}: {e}')
                    handler_results.append((handler_name, 'error', str(e)))
            
            # Log handler results for debugging
            if handler_results:
                results_str = ', '.join([f'{name}:{status}' for name, status, _ in handler_results])
                self.log(f'Error handler results: {results_str}', n=3)
            
            if should_resubmit:
                self.resubmit()
            elif handler_results:
                # All handlers were tried but none recommended retry
                self.log('All error handlers processed, none recommended retry', n=2)
    
    def progress(self, dependency_id=None):
        """Override to handle conditionals, loops, and merge strategies."""
        # Check execution-time conditional
        if self.conditional_enabled and not self.skipped:
            if hasattr(self.__class__, '_conditional_execution_time'):
                if self.__class__._conditional_execution_time:  # Execution-time
                    if not self.evaluate_condition():
                        self.skip()
                        return
        
        # Handle merge strategies for multiple dependencies
        # This must happen BEFORE calling super().progress() because
        # super() checks if wait_ids is empty
        if dependency_id is not None and len(self.dependency_ids) > 1:
            if not hasattr(self, '_completed_dependencies'):
                self._completed_dependencies = set()
            self._completed_dependencies.add(dependency_id)
            
            if self.merge_strategy == 'first':
                # First dependency completed - clear all waits and proceed
                # Remove the completed one first (parent will do this, but we clear all)
                if dependency_id in self.wait_ids:
                    self.wait_ids.remove(dependency_id)
                self.wait_ids.clear()  # Clear all to proceed immediately
            elif self.merge_strategy == 'all':
                # For 'all' strategy, wait for all dependencies
                # If selector is provided, use it to select which dependency's results to use
                # Remove completed dependency from wait_ids
                if dependency_id in self.wait_ids:
                    self.wait_ids.remove(dependency_id)
                # Check if all completed
                if len(self._completed_dependencies) == len(self.dependency_ids):
                    # All completed - if selector is set, use it to choose which dependency to use
                    if self.merge_selector is not None:
                        # Get all completed simulations
                        completed_sims = []
                        for dep_id in self._completed_dependencies:
                            if dep_id in self.dependencies:
                                dep_sim = self.dependencies[dep_id].sim
                                # Selector will filter for completed status
                                completed_sims.append(dep_sim)
                        
                        if completed_sims:
                            try:
                                # Selector should filter for completed and return one
                                selected = self.merge_selector(completed_sims)
                                # Store the selected dependency ID for use in get_dependencies()
                                if selected is not None:
                                    # Find the dependency ID for the selected simulation
                                    for dep_id, dep in self.dependencies.items():
                                        if dep.sim is selected:
                                            self._selected_dependency_id = dep_id
                                            break
                                # Use selected simulation's results when get_dependencies() is called
                                self.wait_ids.clear()
                            except Exception as e:
                                self.warn(f'Error in merge_selector: {e}')
                                self.wait_ids.clear()
                    else:
                        # No selector - standard AND behavior, all dependencies are used
                        # wait_ids will be cleared by parent progress() when empty
                        pass
        
        # Handle loops
        if self.loop_enabled:
            if self.loop_condition is not None:
                if not self.loop_condition(self):
                    # Loop condition false, exit loop
                    self.loop_enabled = False
                elif self.loop_max_iterations is not None:
                    if self.iteration_count >= self.loop_max_iterations:
                        self.loop_enabled = False
        
        # Call parent progress (handles standard AND logic and dependency removal)
        super().progress(dependency_id)
        
        # Handle loop iteration
        if self.loop_enabled and self.finished and not self.failed:
            if self.loop_condition is None or self.loop_condition(self):
                self.next_iteration()
    
    def next_iteration(self):
        """
        Prepare simulation for next loop iteration.
        
        This method:
        1. Increments iteration_count
        2. Resets indicators for next iteration
        3. Calls loop_modifier (if provided) to modify input/state
        """
        self.iteration_count += 1
        
        # Reset indicators for next iteration
        self.reset_indicators()
        
        # Call loop_modifier if provided (allows modifying input at each iteration)
        if self.loop_modifier is not None:
            try:
                self.loop_modifier(self)
            except Exception as e:
                self.warn(f'Error in loop_modifier: {e}')
        
        # Update loop variables if needed
        # (subclasses can override to modify state)
    
    def depends(self, *dependencies, strategy = 'all', selector=None):
        """
        Override to add dependencies with merge strategy support.
        
        Args:
            *dependencies: Dependency tuples (sim, result_name, ...)
            strategy: Merge strategy for multiple dependencies:
                - 'all' (default): Wait for all dependencies (AND logic)
                - 'first': Run when first dependency completes (race condition)
            selector: 
                This function should take a list of completed Simulation objects and return the selected one.
                Example: 
                    selector=lambda sims: min(sims, key=lambda s: s.energy)
        """
        # Store selector for merge strategy
        self.merge_strategy = strategy
        if selector is not None:
            self.merge_selector = selector
        else:
            # Default selector: return first completed simulation
            # This ensures we only work with simulations that have finished
            def completed_sims_selector(sims):
                """Default selector: returns first completed simulation."""
                if not sims:
                    return None
                # Filter to only completed simulations
                completed_sims = [s for s in sims if hasattr(s, 'finished')]
                # If none completed yet, return first one (will be checked later)
                return completed_sims
            self.merge_selector = completed_sims_selector
        
        # Call parent method to add dependencies
        super().depends(*dependencies)
    
    def get_dependencies(self):
        """
        Override to support merge_selector by filtering dependencies.
        
        When a merge_selector has selected a specific dependency,
        only that dependency's results are used.
        """
        # If a dependency was selected via merge_selector, filter to only that one
        if self._selected_dependency_id is not None and self._selected_dependency_id in self.dependencies:
            # Temporarily store original dependencies
            original_dependencies = self.dependencies
            # Create filtered dependencies dict with only the selected one
            from developer import obj
            filtered_deps = obj()
            filtered_deps[self._selected_dependency_id] = self.dependencies[self._selected_dependency_id]
            # Temporarily replace dependencies
            self.dependencies = filtered_deps
            try:
                # Call parent method with filtered dependencies
                super().get_dependencies()
            finally:
                # Restore original dependencies
                self.dependencies = original_dependencies
        else:
            # No selector or no selection made - use all dependencies (standard behavior)
            super().get_dependencies()


def make_enhanced(sim: Simulation, **enhanced_kwargs) -> EnhancedSimulation:
    """
    Convert any Simulation instance to an EnhancedSimulation.
    
    This is a fully generic function that works with any Simulation subclass
    (Pwscf, Qmcpack, Vasp, etc.) without needing type-specific implementations.
    
    Args:
        sim: Any Simulation instance to convert
        **enhanced_kwargs: Enhanced workflow arguments:
            - error_handlers: List of error handlers
            - condition: Conditional function
            - required_machine: Required machine name
            - allowed_machines: Set of allowed machine names
            - loop_enabled: Enable loop support
            - loop_condition: Loop condition function
            - loop_max_iterations: Maximum loop iterations
            - loop_modifier: Function to modify input/state before each iteration
    
    Returns:
        A new instance that inherits from both EnhancedSimulation and sim's class
    
    Example:
        scf = generate_pwscf(...)
        scf_enhanced = make_enhanced(scf, error_handlers=[...], loop_max_iterations=5)
    """
    # Check if already an EnhancedSimulation
    if isinstance(sim, EnhancedSimulation):
        # Just update enhanced kwargs
        for key, value in enhanced_kwargs.items():
            if hasattr(sim, key):
                setattr(sim, key, value)
        return sim
    
    # Create dynamic class that inherits from both EnhancedSimulation and original class
    class_name = f"Enhanced{sim.__class__.__name__}"
    
    # Check if class already exists in module (avoid creating duplicates)
    import sys
    module = sys.modules[__name__]
    if hasattr(module, class_name):
        EnhancedClass = getattr(module, class_name)
    else:
        EnhancedClass = type(class_name, (EnhancedSimulation, sim.__class__), {})
        # Store in module for reuse
        setattr(module, class_name, EnhancedClass)
    
    # Create new instance without calling __init__ (to avoid re-initialization)
    enhanced_sim = EnhancedClass.__new__(EnhancedClass)
    
    # Copy all attributes from original simulation
    # Save fake_sim status first (will be checked/updated after)
    enhanced_sim.__dict__.update(sim.__dict__)
    
    # Initialize ALL enhanced attributes with defaults
    # This ensures all attributes exist even if not in enhanced_kwargs
    # Note: fake_sim is preserved from original (will be updated if machine incompatible)
    enhanced_sim.attempt_number = 0
    enhanced_sim.base_identifier = enhanced_sim.identifier
    enhanced_sim.base_path = enhanced_sim.path
    enhanced_sim.iteration_count = 0
    enhanced_sim.error_handlers = []
    enhanced_sim.attempt_history = []
    enhanced_sim.skipped = False
    enhanced_sim.loop_variables = {}
    enhanced_sim.condition = None
    enhanced_sim.conditional_enabled = False
    enhanced_sim.required_machine = None
    enhanced_sim.allowed_machines = set()
    enhanced_sim.loop_enabled = False
    enhanced_sim.loop_condition = None
    enhanced_sim.loop_max_iterations = None
    enhanced_sim.loop_modifier = None
    enhanced_sim.merge_strategy = 'all'
    enhanced_sim.merge_selector = None
    enhanced_sim._completed_dependencies = set()
    enhanced_sim._selected_dependency_id = None
    
    # Apply enhanced kwargs (overrides defaults)
    enhanced_keys = [
        'error_handlers', 'condition', 
        'required_machine', 'allowed_machines', 'loop_enabled',
        'loop_condition', 'loop_max_iterations', 'loop_modifier',
        'merge_strategy', 'merge_selector'
    ]
    
    for key in enhanced_keys:
        if key in enhanced_kwargs:
            value = enhanced_kwargs[key]
            setattr(enhanced_sim, key, value)
            # Special handling for allowed_machines (convert to set if list)
            if key == 'allowed_machines' and isinstance(value, (list, tuple)):
                enhanced_sim.allowed_machines = set(value)
    
    # Set conditional_enabled based on condition (after applying kwargs)
    if enhanced_sim.condition is not None:
        enhanced_sim.conditional_enabled = True
    
    # Helper function to replace original in all_sims and update dependency references
    def _replace_in_all_sims():
        if hasattr(Simulation, 'all_sims'):
            # Find and replace the original simulation in all_sims
            # Use identity comparison (is) to avoid recursion issues with __eq__
            for i, existing_sim in enumerate(Simulation.all_sims):
                if existing_sim is sim:
                    Simulation.all_sims[i] = enhanced_sim
                    break
            
            # Also update any dependency references in other simulations
            # This is critical: if scf3 depends on scf2_base, we need to update it to scf2
            for other_sim in Simulation.all_sims:
                if other_sim is not enhanced_sim and hasattr(other_sim, 'dependencies'):
                    # Check dependencies dict
                    for dep_id, dep in list(other_sim.dependencies.items()):
                        if hasattr(dep, 'sim') and dep.sim is sim:
                            dep.sim = enhanced_sim
                            # Update the key if simid changed (though it shouldn't)
                            if dep_id != enhanced_sim.simid:
                                other_sim.dependencies[enhanced_sim.simid] = dep
                                del other_sim.dependencies[dep_id]
                                other_sim.dependency_ids.discard(dep_id)
                                other_sim.dependency_ids.add(enhanced_sim.simid)
                                other_sim.wait_ids.discard(dep_id)
                                other_sim.wait_ids.add(enhanced_sim.simid)
                    
                    # Check dependents dict
                    if hasattr(other_sim, 'dependents'):
                        for dep_id, dependent in list(other_sim.dependents.items()):
                            if dependent is sim:
                                other_sim.dependents[enhanced_sim.simid] = enhanced_sim
                                if dep_id != enhanced_sim.simid:
                                    del other_sim.dependents[dep_id]
                    
                    # Update dependency_ids and wait_ids sets
                    if hasattr(other_sim, 'dependency_ids') and sim.simid in other_sim.dependency_ids:
                        other_sim.dependency_ids.discard(sim.simid)
                        other_sim.dependency_ids.add(enhanced_sim.simid)
                    if hasattr(other_sim, 'wait_ids') and sim.simid in other_sim.wait_ids:
                        other_sim.wait_ids.discard(sim.simid)
                        other_sim.wait_ids.add(enhanced_sim.simid)
    
    # Check machine compatibility IMMEDIATELY after setting required_machine/allowed_machines
    # This is critical since we bypassed __init__ where this normally happens
    if enhanced_sim.required_machine or enhanced_sim.allowed_machines:
        # Get current machine for explicit comparison
        machine_name = None
        if hasattr(enhanced_sim, 'job') and enhanced_sim.job is not None:
            machine_name = getattr(enhanced_sim.job, 'machine', None)
        if machine_name is None:
            from machines import Job
            machine_name = Job.machine
        if machine_name is None:
            from project_manager import ProjectManager
            if ProjectManager.machine is not None:
                machine_name = getattr(ProjectManager.machine, 'name', None)
        
        # Explicit check: if required_machine is set, verify it matches
        if enhanced_sim.required_machine:
            if machine_name is None:
                # Can't determine machine - mark as fake to be safe
                enhanced_sim.fake_sim = True
                _replace_in_all_sims()
                return enhanced_sim
            elif str(machine_name).strip().lower() != str(enhanced_sim.required_machine).strip().lower():
                # Machine doesn't match - mark as fake
                enhanced_sim.fake_sim = True
                _replace_in_all_sims()
                return enhanced_sim
        
        # Explicit check: if allowed_machines is set, verify machine is in list
        if enhanced_sim.allowed_machines and len(enhanced_sim.allowed_machines) > 0:
            if machine_name is None:
                # Can't determine machine - mark as fake to be safe
                enhanced_sim.fake_sim = True
                _replace_in_all_sims()
                return enhanced_sim
            else:
                machine_normalized = str(machine_name).strip().lower()
                allowed_normalized = {str(a).strip().lower() for a in enhanced_sim.allowed_machines}
                if machine_normalized not in allowed_normalized:
                    # Machine not in allowed list - mark as fake
                    enhanced_sim.fake_sim = True
                    _replace_in_all_sims()
                    return enhanced_sim
        
        # Also call is_machine_compatible for consistency (it will also set fake_sim)
        compatible = enhanced_sim.is_machine_compatible()
        if not compatible:
            enhanced_sim.fake_sim = True
            _replace_in_all_sims()
            return enhanced_sim
    
    # Definition-time conditional check (if applicable)
    if enhanced_sim.conditional_enabled:
        cls = enhanced_sim.__class__
        if hasattr(cls, '_conditional_execution_time'):
            if not cls._conditional_execution_time:  # Definition-time
                if not enhanced_sim.evaluate_condition():
                    enhanced_sim.fake_sim = True
                    _replace_in_all_sims()
                    return enhanced_sim
    
    # CRITICAL: Replace original simulation in Simulation.all_sims with enhanced version
    # This ensures that when add_simulations() is called without arguments,
    # it uses the enhanced version (which may be blocked) instead of the original
    _replace_in_all_sims()
    
    return enhanced_sim


def create_branch(parent: Simulation, branches: List[Tuple[Simulation, Optional[Callable], str]], **enhanced_kwargs) -> List[EnhancedSimulation]:
    """
    Create conditional branches from a parent simulation.
    
    This provides a natural way to create if/else or multi-way conditional branches
    without needing separate make_enhanced() calls for each branch.
    
    Args:
        parent: Parent simulation that branches depend on
        branches: List of (sim_base, condition, result_name) tuples where:
            - sim_base: Base simulation instance (before make_enhanced)
            - condition: Optional conditional function (Callable[[Simulation], bool])
                         If None, branch always runs (unconditional)
            - result_name: Result name from parent (e.g., 'charge_density')
        **enhanced_kwargs: Additional enhanced kwargs applied to all branches
                          (error_handlers, required_machine, etc.)
    
    Returns:
        List of EnhancedSimulation instances (one for each branch)
    
    Example:
        # Create if/else branch
        scf3, scf4 = create_branch(
            parent=scf2,
            branches=[
                (scf3_base, energy_below_threshold, 'charge_density'),
                (scf4_base, energy_above_threshold, 'charge_density'),
            ],
            error_handlers=[RetryErrorHandler(max_retries=2)]
        )
        
        # Create multi-way branch
        scf_a, scf_b, scf_c = create_branch(
            parent=scf1,
            branches=[
                (scf_a_base, condition_a, 'charge_density'),
                (scf_b_base, condition_b, 'charge_density'),
                (scf_c_base, None, 'charge_density'),  # Always runs
            ]
        )
    """
    enhanced_branches = []
    
    for i, branch_spec in enumerate(branches):
        if len(branch_spec) < 3:
            raise ValueError(f'Branch {i} must be (sim_base, condition, result_name) tuple')
        
        sim_base, condition, result_name = branch_spec[0], branch_spec[1], branch_spec[2]
        
        # Ensure sim_base depends on parent
        if not hasattr(sim_base, 'dependencies') or parent.simid not in sim_base.dependency_ids:
            sim_base.depends((parent, result_name))
        
        # Create enhanced simulation with condition
        branch_kwargs = enhanced_kwargs.copy()
        if condition is not None:
            branch_kwargs['condition'] = condition
        
        enhanced_branch = make_enhanced(sim_base, **branch_kwargs)
        enhanced_branches.append(enhanced_branch)
    
    return enhanced_branches

