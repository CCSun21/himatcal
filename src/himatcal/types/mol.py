from __future__ import annotations

from pydantic import BaseModel, field_validator


class CASNumber(BaseModel):
    """CAS number validation class (AI-generated).

    Validates Chemical Abstracts Service (CAS) registry numbers according to standard format:
    - 2-7 digits followed by a hyphen
    - Then 2 digits and another hyphen
    - Finally 1 digit
    - Total digits cannot exceed 10
    """

    cas_number: str

    @field_validator("cas_number")
    def validate_cas_number(cls, value):
        import re

        """Validate CAS number format and length."""
        # Validate basic format
        pattern = re.compile(r"^\d{2,7}-\d{2}-\d{1}$")
        if not re.match(pattern, value):
            raise ValueError(
                "Invalid CAS number format. It should be 2-7 digits followed by a hyphen, then 2 digits, and another hyphen followed by 1 digit."
            )

        # Validate first part length
        parts = value.split("-")
        if len(parts[0]) > 7:  # Check if first part exceeds 7 digits
            raise ValueError("Invalid CAS number: first part cannot exceed 7 digits")

        # Validate total length
        total_digits = sum(len(part) for part in parts)
        if total_digits > 10:  # Total digits in CAS number should not exceed 10
            raise ValueError("Invalid CAS number: total length cannot exceed 10 digits")

        return value
