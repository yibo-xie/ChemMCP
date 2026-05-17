import argparse
import functools
import inspect
import logging
from typing import List, Optional, Dict
from mcp.server.fastmcp import FastMCP

from .errors import catch_errors

logger = logging.getLogger(__name__)

mcp_instance = FastMCP("ChemMCP")


original_tool = mcp_instance.tool

def tool_with_catch(*dargs, **dkwargs):
    def decorator(fn):
        # first wrap errors, then register as a tool
        wrapped = catch_errors(fn)
        return original_tool(*dargs, **dkwargs)(wrapped)
    return decorator

mcp_instance.tool = tool_with_catch


class ChemMCPManager:
    _tools: Dict[str, type] = {}

    @staticmethod
    def register_tool(cls):
        cls._registered_mcp_tool = True
        ChemMCPManager._tools[cls.name] = cls
        return cls
    
    @staticmethod
    def get_registered_tools():
        return list(ChemMCPManager._tools.keys())
    
    @staticmethod
    def discover_and_register():
        """Import all tool modules to trigger @register_tool decorators.
        
        The tools package uses lazy imports (__getattr__), so tool classes
        are not registered until their modules are actually imported. This
        method walks the public attributes of chemmcp.tools to force-load
        every tool module, ensuring all @register_tool decorators execute.
        
        Returns:
            int: Number of tools discovered and registered.
        """
        import chemmcp.tools
        count_before = len(ChemMCPManager._tools)
        for name in dir(chemmcp.tools):
            if not name.startswith('_'):
                try:
                    getattr(chemmcp.tools, name)
                except Exception:
                    logger.debug(f"Skipping tool {name} due to import error")
        count_after = len(ChemMCPManager._tools)
        discovered = count_after - count_before
        logger.info(f"Discovered {discovered} tools (total: {count_after})")
        return count_after
    
    @staticmethod
    def init_mcp_tools(mcp, tools: Optional[List[str]] = None):
        registered_tool_names = set(ChemMCPManager.get_registered_tools())
        
        if tools is None:
            tools = registered_tool_names
        else:
            tools = set(tools)
            for tool_name in tools:
                if tool_name not in registered_tool_names:
                    raise ValueError(f"Tool '{tool_name}' does not exist.")
                
        num_tools = len(tools)
        initialized = 0
        skipped = []
        for tool_name in tools:
            cls = ChemMCPManager._tools[tool_name]

            try:
                inst = cls()

                sig = inspect.signature(cls._run_base)
                params = list(sig.parameters.values())[1:]
                exposed_sig = sig.replace(parameters=params)

                def make_wrapper(inst, func):
                    @functools.wraps(func)
                    def wrapper(*args, **kwargs):
                        return inst.run_code(*args, **kwargs)
                    return wrapper

                wrapper = make_wrapper(inst, cls._run_base)
                wrapper.__signature__ = exposed_sig
                wrapper.__name__      = cls.func_name
                wrapper.__doc__       = cls.get_doc(interface='code')

                existing = set(mcp._tool_manager._tools.keys())
                if cls.func_name not in existing:
                    mcp.tool()(wrapper)
                initialized += 1
            except Exception as e:
                logger.warning("Skipping tool '%s': %s", tool_name, e)
                skipped.append(tool_name)
                continue
        
        logger.info(f"Initialized {initialized}/{num_tools} tools to MCP." + (f" Skipped: {skipped}" if skipped else ""))


def run_mcp_server():
    parser = argparse.ArgumentParser(description="Run the MCP server.")
    parser.add_argument('--tools', default=None, type=str, nargs='+', help="The tools to load.")
    parser.add_argument('--sse', action='store_true', help="Run the server with SSE (Server-Sent Events) support.")
    args = parser.parse_args()

    ChemMCPManager.discover_and_register()
    ChemMCPManager.init_mcp_tools(mcp_instance, args.tools)

    logger.info("Initialized ChemMCP tools.")

    if args.sse:
        # build a Starlette/uvicorn app
        app = mcp_instance.sse_app()
        import uvicorn
        uvicorn.run(app, host="127.0.0.1", port=8001)
    else:
        # Run the MCP server with standard input/output
        mcp_instance.run(transport='stdio')
