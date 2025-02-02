defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction
  """
  
  alias Wavelets.DWT
  alias Wavelets.Filter
  alias Wavelets.IDWT

  @doc """
  Performs 1D wavelet packet decomposition
  """
  @spec decompose_1d(list(number), Filter.t(), non_neg_integer(), Wavelets.precision()) ::
    map()
  def decompose_1d(signal, filter, levels, _precision \\ :double) do
    # Start with root node
    tree = %{{0, 0} => signal}
    
    # Build each level
    Enum.reduce(1..levels, tree, fn level, acc ->
      build_level(acc, filter, level)
    end)
  end

  defp build_level(tree, filter, level) do
    prev_level = level - 1
    
    # Get all nodes from previous level
    prev_nodes = for {key = {^prev_level, _}, signal} <- tree, into: %{} do
      {key, signal}
    end
    
    # Create new nodes
    new_nodes = 
      Enum.flat_map(Map.to_list(prev_nodes), fn {{_, j}, signal} ->
        {approx, details} = DWT.forward_1d(signal, filter)
        [
          {{level, j * 2}, approx},
          {{level, j * 2 + 1}, details}
        ]
      end)
      |> Map.new()
    
    Map.merge(tree, new_nodes)
  end

  @doc """
  Reconstructs signal from wavelet packet decomposition
  """
  @spec reconstruct_1d(map(), Filter.t(), Wavelets.precision()) :: list(number)
  def reconstruct_1d(tree, filter, _precision \\ :double) do
    max_level = tree
                |> Map.keys()
                |> Enum.map(fn {level, _} -> level end)
                |> Enum.max()
    
    reconstruct_level(tree, filter, max_level)
  end

  defp reconstruct_level(tree, _filter, 0) do
    Map.get(tree, {0, 0})
  end
  
  defp reconstruct_level(tree, filter, level) do
    nodes = for j <- 0..(:math.pow(2, level-1) |> trunc() |> Kernel.-(1)),
                approx = Map.get(tree, {level, j * 2}),
                details = Map.get(tree, {level, j * 2 + 1}),
                approx != nil and details != nil do
      {{level - 1, j}, IDWT.inverse_1d(approx, details, filter)}
    end
    |> Map.new()

    Map.merge(tree, nodes)
    |> reconstruct_level(filter, level - 1)
  end
end