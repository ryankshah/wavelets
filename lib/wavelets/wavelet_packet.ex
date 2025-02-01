defmodule Wavelets.WaveletPacket do
  @moduledoc """
  Implementation of Wavelet Packet decomposition and reconstruction
  """

  alias Wavelets.{DWT, Filter, IDWT}

  @doc """
  Performs 1D wavelet packet decomposition
  """
  @spec decompose_1d(
          list(number),
          Filter.t(),
          non_neg_integer(),
          Wavelets.precision()
        ) ::
          map()
  def decompose_1d(signal, filter, levels, _precision \\ :double) do
    tree = %{{0, 0} => signal}
    decompose_level(tree, filter, levels)
  end

  defp decompose_level(tree, _filter, 0), do: tree

  defp decompose_level(tree, filter, level) do
    current_level = level - 1
    nodes = get_level_nodes(tree, current_level)

    new_nodes =
      nodes
      |> Enum.flat_map(fn {node_idx, signal} ->
        {approx, details} = DWT.forward_1d(signal, filter)

        [
          {{level, node_idx * 2}, approx},
          {{level, node_idx * 2 + 1}, details}
        ]
      end)
      |> Map.new()

    tree
    |> Map.merge(new_nodes)
    |> decompose_level(filter, level - 1)
  end

  defp get_level_nodes(tree, level) do
    tree
    |> Map.filter(fn {{node_level, _}, _} -> node_level == level end)
    |> Map.filter(fn {_, signal} -> is_list(signal) end)
  end

  @doc """
  Reconstructs signal from wavelet packet decomposition
  """
  @spec reconstruct_1d(map(), Filter.t(), Wavelets.precision()) :: list(number)
  def reconstruct_1d(tree, filter, _precision \\ :double) do
    max_level =
      tree
      |> Map.keys()
      |> Enum.map(fn {level, _} -> level end)
      |> Enum.max()

    reconstruct_level(tree, filter, max_level)
  end

  defp reconstruct_level(tree, _filter, 0) do
    Map.get(tree, {0, 0})
  end

  defp reconstruct_level(tree, filter, level) do
    max_nodes = trunc(:math.pow(2, level - 1))

    new_nodes =
      for j <- 0..(max_nodes - 1),
          approx = Map.get(tree, {level, j * 2}),
          details = Map.get(tree, {level, j * 2 + 1}),
          approx != nil and details != nil do
        {{level - 1, j}, IDWT.inverse_1d(approx, details, filter)}
      end
      |> Map.new()

    tree
    |> Map.merge(new_nodes)
    |> reconstruct_level(filter, level - 1)
  end
end
