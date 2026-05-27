# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Haohao Fu <fhh2626@nankai.edu.cn, fhh2626@gmail.com>

import pathlib
from PySide6.QtGui import QIcon

def get_app_icon():
    """Get the application icon.
    
    Returns:
        QIcon: The application icon, or a default icon if the file is not found.
    """
    icon_path = pathlib.Path(__file__).parent / "resources" / "icon.png"
    if icon_path.exists():
        return QIcon(str(icon_path))
    return QIcon()  # Return empty icon if not found
