from datetime import datetime
import logging
from pathlib import Path


class ColoredFormatter(logging.Formatter):
    COLORS = {
        logging.DEBUG: "\033[34m",  # BLUE
        logging.INFO: "\033[32m",  # GREEN
        logging.WARNING: "\033[38;5;214m",  # ORANGE
        logging.ERROR: "\033[31m",  # RED
        logging.CRITICAL: "\033[1;31m",  # BOLD RED
    }
    MODULE_COLOR = "\033[38;5;45m"  # Neon blue

    def format(self, record):
        RESET = "\033[0m"
        color = self.COLORS.get(record.levelno, RESET)
        record.levelname = f"{color}{record.levelname}{RESET}"
        record.module = f"{self.MODULE_COLOR}{record.module}{RESET}"
        format_string = "| %(levelname)s | %(module)s: - %(message)s"
        formatter = logging.Formatter(format_string)
        return formatter.format(record)


class Logger:
    _logger_instance = None 
    _log_path = None

    @classmethod
    def setup_logger(cls, log_path="sip.log"):
        logger = logging.getLogger("sip_logger")
        logger.setLevel(logging.INFO)

        timestamp = datetime.now().strftime("%d%m%y_%H%M%S")

        original_path = Path(log_path)
        new_stem = f"{original_path.stem}_{timestamp}"
        log_path = original_path.with_name(new_stem).with_suffix(original_path.suffix)

        file_handler = logging.FileHandler(log_path, encoding='utf-8')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(
            logging.Formatter(
                "%(asctime)s | %(levelname)s | %(module)s: - %(message)s"
            )
        )

        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(ColoredFormatter())

        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

        cls._logger_instance = logger
        cls._log_path = log_path

        return log_path
    

    @classmethod
    def get_logger(cls, log_path="sip.log"): #, clear_logs=False, queue=None):
        if cls._logger_instance is None:
            log_path = cls.setup_logger(log_path=log_path)#, clear_logs=clear_logs, queue=queue)
        return cls._logger_instance
