from __future__ import annotations

import re

from pydantic import BaseModel, field_validator


class CASNumber(BaseModel):
    cas_number: str

    @field_validator("cas_number")
    def validate_cas_number(cls, value):
        """Validate CAS number format and length."""
        # 验证基本格式
        pattern = re.compile(r"^\d{2,7}-\d{2}-\d{1}$")
        if not re.match(pattern, value):
            raise ValueError(
                "Invalid CAS number format. It should be 2-7 digits followed by a hyphen, then 2 digits, and another hyphen followed by 1 digit."
            )

        # 验证第一部分的长度
        parts = value.split("-")
        if len(parts[0]) > 7:  # 检查第一部分是否超过7位数字
            raise ValueError("Invalid CAS number: first part cannot exceed 7 digits")

        # 验证整体长度
        total_digits = sum(len(part) for part in parts)
        if total_digits > 10:  # CAS号的数字总数不应超过10位
            raise ValueError("Invalid CAS number: total length cannot exceed 10 digits")

        return value
