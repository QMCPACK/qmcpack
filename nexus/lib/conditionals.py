"""
Conditional utilities for enhanced workflow simulations.

Provides utility functions for conditional execution and machine filtering.
"""

from typing import Callable, Optional, Any, Dict


def _add_conditional_metadata(func: Callable, description: str, metadata: Optional[Dict[str, Any]] = None) -> Callable:
    """
    Add metadata to a conditional function for better status reporting.
    
    Args:
        func: Conditional function to annotate
        description: Human-readable description of the conditional
        metadata: Additional metadata dictionary
        
    Returns:
        Function with metadata attached
    """
    func._conditional_description = description
    func._conditional_metadata = metadata or {}
    return func


def get_conditional_description(condition: Optional[Callable]) -> str:
    """
    Get a human-readable description of a conditional function.
    
    Args:
        condition: Conditional function (may have metadata attached)
        
    Returns:
        Description string, or default if no metadata available
    """
    if condition is None:
        return "None"
    
    # Handle non-callable objects gracefully
    if not callable(condition):
        return f"<non-callable:{type(condition).__name__}>"
    
    try:
        # Check for attached metadata
        if hasattr(condition, '_conditional_description'):
            return condition._conditional_description
        
        # Try to infer from function name
        if hasattr(condition, '__name__'):
            name = condition.__name__
            if name != '<lambda>':
                return name
        
        # For lambda functions, try to get source if available
        if hasattr(condition, '__code__'):
            # Could extract more info from __code__ if needed
            pass
        
        # Default fallback
        return "custom"
    except Exception:
        # If anything goes wrong, return a safe fallback
        return "custom"


def machine_conditional(machine_name: str, *, required: bool = True) -> Callable:
    """
    Create a conditional function that checks machine compatibility.
    
    Args:
        machine_name: Name of machine to check
        required: If True, must match; if False, must not match
        
    Returns:
        Conditional function
    """
    def check_machine(sim: Any) -> bool:
        current_machine = getattr(sim, 'job', None)
        if current_machine is None:
            return not required
        
        machine_name_attr = getattr(current_machine, 'machine', None)
        if machine_name_attr is None:
            return not required
        
        matches = machine_name_attr == machine_name
        return matches if required else not matches
    
    description = f"machine={'==' if required else '!='}{machine_name}"
    return _add_conditional_metadata(check_machine, description, {'machine_name': machine_name, 'required': required})


def result_conditional(dependency: str, check_func: Callable[[Any], bool]) -> Callable:
    """
    Create a conditional function that checks dependency results.
    
    Args:
        dependency: Name of dependency to check
        check_func: Function that takes dependency result and returns bool
        
    Returns:
        Conditional function
    """
    def check_result(sim: Any) -> bool:
        if not hasattr(sim, 'dependencies'):
            return False
        
        for dep in sim.dependencies.values():
            dep_sim = dep.sim
            if hasattr(dep_sim, 'outputs') and dep_sim.outputs is not None:
                result = getattr(dep_sim.outputs, dependency, None)
                if result is not None:
                    return check_func(result)
        
        return False
    
    check_func_name = getattr(check_func, '__name__', 'custom')
    description = f"result:{dependency}({check_func_name})"
    return _add_conditional_metadata(check_result, description, {'dependency': dependency, 'check_func': check_func_name})


def combine_conditionals(*conditionals: Callable, 
                        operator: str = 'and') -> Callable:
    """
    Combine multiple conditionals with AND or OR logic.
    
    Args:
        *conditionals: Conditional functions to combine
        operator: 'and' or 'or'
        
    Returns:
        Combined conditional function
    """
    allowed_operators = ['and', 'or']
    if operator not in allowed_operators:
        raise ValueError(f"Unknown operator: {operator}. Use one of: {allowed_operators}")
    
    def combined(sim: Any) -> bool:
        if operator == 'and':
            return all(cond(sim) for cond in conditionals)
        elif operator == 'or':
            return any(cond(sim) for cond in conditionals)
    
    cond_descriptions = [get_conditional_description(cond) for cond in conditionals]
    op_str = f' {operator.upper()} '
    description = f"({op_str.join(cond_descriptions)})"
    return _add_conditional_metadata(combined, description, {'operator': operator, 'conditionals': conditionals})


