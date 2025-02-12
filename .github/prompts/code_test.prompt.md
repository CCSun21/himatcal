# Test File Generation Prompt

When I need to create a test file for a Python file, please follow these guidelines:

## File Naming Convention

1. The test file should be named test_<original_filename>.py
2. The test file should be placed in the tests folder in the project's root directory

## Testing Framework and Tools

1. Use the following tools and frameworks:
2. pytest as the primary testing framework
3. pytest-cov for code coverage reports
4. unittest.mock for mocking external dependencies

## Test Structure Requirements

1. Each test class should correspond to a class being tested
2. The test class should be named Test<OriginalClassName>

## Include the following types of tests:

1. Unit tests: Testing individual functions/methods
2. Parameterized tests: Testing different input scenarios
3. Exception tests: Validating error handling
4. Boundary tests: Validating boundary conditions

## Test Content Requirements

Each test should:

1. Include a clear docstring describing the test's purpose
2. Follow the Arrange-Act-Assert pattern
3. Use meaningful test data
4. Validate expected output and side effects
5. Handle all possible err

## Code example

```python
import pytest
from unittest.mock import Mock, patch
from your_module import YourClass

class TestYourClass:
    """Test suite for YourClass."""
    
    def setup_method(self):
        """Set up test fixtures before each test method."""
        self.instance = YourClass()
    
    def test_method_name(self):
        """Test specific functionality."""
        # Arrange
        input_data = "test"
        expected = "result"
        
        # Act
        result = self.instance.method(input_data)
        
        # Assert
        assert result == expected
    
    @pytest.mark.parametrize("input,expected", [
        ("case1", "result1"),
        ("case2", "result2"),
    ])
    def test_method_with_different_inputs(self, input, expected):
        """Test method with different input values."""
        assert self.instance.method(input) == expected
    
    def test_method_raises_exception(self):
        """Test exception handling."""
        with pytest.raises(ValueError):
            self.instance.method(invalid_input)
```

## Guidelines for Test Writing

1. Ensure test independence
2. Avoid dependencies between tests
3. Use fixtures and setup/teardown judiciously
4. Keep tests simple and clear
5. Isolate external dependencies with mock objects when appropriate

## Test Coverage Targets

1. Line coverage: At least 80%
2. Branch coverage: At least 80%
3. Path coverage: Cover as many critical paths as possible

## Execution Commands

```bash
Run tests
pytest tests/

Generate coverage report
pytest --cov=your_package tests/
```
